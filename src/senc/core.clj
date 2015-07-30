(ns senc.core
  (:require [clojure.core.matrix :as mx]
            [schema.core :as s]
            [senc.estimate.mle :as mle]
            [clojure.java.io :as io]
            [clojure.set :as set]
            [clojure.core.matrix.impl.pprint :as mpp :refer [pm]]
            [clojure.pprint :refer [pprint]]
            [munge.core])
  (:use [clojure.tools.cli :only [cli]]
        [senc.schema :only [Vec Mat ProbVec BinVec BinMat]]
        [munge.io.matrix-mm :only [load-matrix save-matrix]]
        [munge.io.data-frame :only [save-data-frame]]
        [munge.io.core :only [load-ids]]
        [munge.matrix :only [proportional select-rows selected-rows
                             binary-vector b-or vec-non-zeros
                             create-sparse-indexed-vector
                             sparse-indexed-vector
                             sparse-row-matrix
                             sparse-column-matrix]]
        [senc.util :only [safe-log log-multicat]])
  (:import [mikera.matrixx.impl SparseRowMatrix SparseColumnMatrix]
           [mikera.vectorz.impl SparseIndexedVector SparseHashedVector ZeroVector]
           [mikera.vectorz Vectorz]
           [mikera.vectorz.util DoubleArrays]
           [mikera.indexz Index])
  (:gen-class))

(set! *warn-on-reflection* true)
(mx/set-current-implementation :vectorz)

(s/defn get-private-field :- s/Any
  [x :- s/Any
   field :- s/Str]
  (let [^java.lang.reflect.Field f (-> ^java.lang.Object x
                                       (.getClass)
                                       (.getDeclaredField field))]
    (.setAccessible f true)
    (.get f x)))

(defn print-matrix-info
  [msg m]
  (println msg
           (type m)
           (mx/sparse? m)
           (when (or (seq? m) (instance? java.util.List m))
             (format "seq contains: %s" (type (first m)))))
  m)

(defn write-log
  [msg]
  (println (format "%s - %s" (str (java.util.Date.)) msg)))

(s/defn round-to-zero! :- (s/either Mat Vec)
  "Returns a new matrix (need to recompute sparse index) with 
  values below threshold set to 0."
  [threshold :- s/Num
   m :- (s/either Mat Vec)]
  ;; TODO: best speed assumes row matrix, fixable through current core.matrix API?
  (if (mx/vec? m)
    (.roundToZero ^SparseIndexedVector m threshold)
    (let [vs (if (mx/vec? m) [m] (mx/rows m))]
      (SparseRowMatrix/create ^java.util.List (map (fn [v] (when (instance? SparseIndexedVector v)
                                                             (-> (.roundToZero ^SparseIndexedVector v threshold)
                                                                 (proportional))))
                                                   vs)))))

(s/defn restricted-objs :- Mat
  "A convenience function for restricting object feature matrix to
  those objects belonging to the community's seed subgraph."
  [seed-groups-objs :- BinMat
   objs-feats :- Mat
   comm-idx :- s/Int]
  (->> (mx/get-row seed-groups-objs comm-idx)
       (selected-rows objs-feats)))

(s/defn restricted-comms :- Mat
  "A convenience function for restricting the community params
  matrix according to those communities influential to a set of objects."
  [objs-comms :- BinMat
   comms-props :- Mat
   obj-idxs :- #{s/Int}]
  (let [row-indicator (apply b-or (select-rows objs-comms (vec obj-idxs)))]
    (selected-rows comms-props row-indicator)))

(s/defn estimate-props :- Mat
  "Used to calculate the initial estimates for community feature proportions."
  [objs-feat-vals :- Mat
   groups-objs :- BinMat]
  ;; Was using this instead of doseq:
  (let [ncols (mx/column-count objs-feat-vals)
        result (sparse-row-matrix (->> (pmap
                                        (fn [group-objs]
                                          (let [accum (SparseIndexedVector/createLength ncols)]
                                            (doseq [row (mx/rows (selected-rows objs-feat-vals group-objs))]
                                              (when-not (mx/zero-matrix? row)
                                                (mx/add! accum row)))
                                            (-> accum proportional)))
                                        (mx/rows groups-objs))
                                       (into [])))]
    result))

;; TODO: pass in a column-matrix version of comm-sparams for speed up?
;; TODO: round-to-zero! call at end to help performance in later ops that use this output?
(s/defn estimate-memb :- Vec
  "Note: use a filtered comm-params matrix if you want to force zero membership 
  in some communities."
  [obj-feats :- Vec
   comms-params :- Mat]
  ;; Calculate probability of object from weighted counts for all communities

  ;; TODO: we don't need all columns, how can we only do calculations with non-zero columns?
  
  (->> (mx/columns comms-params)
       (map mx/mutable)
       (map proportional)
       (map (s/fn [term-count :- s/Num
                    weight-col :- Vec]                 
              (mx/emul! weight-col term-count)
              weight-col)
            obj-feats)
       ;; TODO: can remove vec? removing vec makes it really slooow. coll -> array faster than seq -> array?
       (vec)
       (sparse-column-matrix)
       ;; TODO: rows from column matrix, expensive!!!
       (mx/rows)
       (map mx/mutable)
       ;; TODO: is there a map we can use that will generate a sparse vector?
       (map mx/esum)
       ;; TODO: can remove vec?
       (vec)
       (sparse-indexed-vector)
       (proportional)))

;; TODO: convert into a function that can be used for computing single object memberships?
(s/defn estimate-membs :- Mat
  "Estimate the membership of a group of objects, specified by index. Used to limit estimation
  to only those objects belonging to seed groups."
  [objs-feat-vals :- Mat
   objs-groups :- BinMat
   comms-props :- Mat
   obj-idxs :- [s/Int]]
  ;; TODO: can we convert comms-props to columns here to avoid
  ;;       conversion in every call to estimate-memb?
  (let [num-objs (mx/row-count objs-feat-vals)
        num-comms (mx/row-count comms-props)
        membs (SparseRowMatrix/create num-objs num-comms)]
    ;; TODO: Fix this below?
    ;; Memory still grows (because of mapping a closure?) but
    ;; doesn't appear to OOM. Fills up heap to max and then does a big cleanup.
    ;; Hypothesis is that now it clean up the mess since we're not holding the
    ;; head. Still annoyed that it gets so big in first place....
    ;; An attempt to not hold the head...
    (doseq [[idx obj-row] (pmap (fn [obj-idx]
                                  [obj-idx (estimate-memb (mx/get-row objs-feat-vals obj-idx)
                                                          (restricted-comms objs-groups
                                                                            comms-props
                                                                            (hash-set obj-idx)))])
                                obj-idxs)]
      (mx/set-row! membs idx obj-row))

    ;; Avoids the memory explosion..
    ;; (dotimes [obj-idx num-objs]
    ;;   (let [obj-row (estimate-memb (mx/get-row objs-feat-vals obj-idx)
    ;;                                 (restricted-comms objs-groups
    ;;                                                   comms-props
    ;;                                                   (hash-set obj-idx)))]
    ;;        ;; TODO: get this non-zeros and row setter under control. core.matrix should support this!!!
    ;;     (mx/set-row! membs obj-idx obj-row)))
    
    membs))

;; TODO: Use the BinMat objs-groups and selected-rows instead of a select-obj-comms fn
;; or make m-step use a similar function. Check which is cleaner.
(s/defn e-step :- Mat
  "Calculate new object membership probabilities.
  The select-comms fn is used to select the communities
  of which an object is a member."
  [select-obj-comms :- (s/=> Mat Mat #{s/Int}) ; This is restricted-comms with obj-comms partially applied
   objs-feats :- Mat
   comms-params :- Mat]
  ;; TODO: Make cleaner? This is a little goofy, building up a seq of comm-param matrices, one for each object by index.
  (let [num-objs (mx/row-count objs-feats)
        obj-comm-params (map (comp (partial select-obj-comms comms-params) hash-set)
                             (range num-objs))
        result (sparse-row-matrix (pmap estimate-memb (mx/rows objs-feats) obj-comm-params))]
    result))

;; phi
(s/defn subgraph-membership :- ProbVec
  "Calculate the community membership weights for the group of objects making up the subgraph."
  [objs-membs :- Mat
   obj-prop-num-obs :- ProbVec]
  (mx/mmul obj-prop-num-obs objs-membs))

;; Theta_c!=i
(s/defn mixed-props :- Vec
  "The subgraph-memb is almost a ProbVec but the element corresponding the community being estimated
  will be zero.

  The result will be not quite a ProbVec either, as the community selected for estimation is not contributing."
  ;; TODO: we could return a proportional vector? But what if all zeros in the trivial case of objects that are
  ;; only members of a single community? Best to treat this result as non-proportional for consistency?
  [subgraph-memb :- Vec
   subset-comm-props :- Mat]
  (mx/mmul subgraph-memb
           subset-comm-props))

(s/defn estimate-comm-props :- ProbVec
  [groups-objs :- BinMat
   objs-groups :- BinMat
   objs-membs :- Mat
   objs-feats :- Mat
   comms-props :- Mat
   comm-idx :- s/Int]
  (let [;;seed-objs (mx/get-row groups-objs comm-idx) ; objs part of the selected community seed subgraph
        num-objs (mx/column-count groups-objs)
        seed-objs (let [objs (mx/get-column objs-membs comm-idx)
                        objs-bin (create-sparse-indexed-vector num-objs)]
                    (doseq [nzi (mx/non-zero-indices objs)]
                      (mx/mset! objs-bin nzi 1))
                    objs-bin)

        ;; comms which objs from selected comm may be part of
        subset-comm-props (selected-rows comms-props
                                         (apply b-or (selected-rows objs-groups seed-objs)))

        ;; TODO: Need to update matrix_api.clj to use SparseIndexedVector.
        ;;       Or it was just that addCopy isn't supported but regular add is fine?
        seed-objs-feats (reduce mx/add! (SparseIndexedVector/createLength (mx/column-count objs-feats))
                                (selected-rows objs-feats seed-objs))

        ;; contribution proportion based on number of features per-object
        ;; TODO: rows and columns conversion here
        obj-prop-num-obs (->> (selected-rows objs-feats seed-objs)
                              mx/columns
                              (reduce mx/add! (SparseIndexedVector/createLength (mx/row-count objs-feats)))
                              proportional)
        ;; the "averaged" phi for the selected community's seed subgraph, this is the complete phi vector
        ;; TODO: should we use round-to-zero to reduce the number of computations performed in mixed-props?
        ;;subgraph-memb (subgraph-membership (selected-rows objs-membs seed-objs) obj-prop-num-obs)
        subgraph-memb (->> (subgraph-membership (selected-rows objs-membs seed-objs) obj-prop-num-obs)
                           (round-to-zero! 1e-4)
                           proportional)
        other-props (mixed-props (mx/mset subgraph-memb comm-idx 0) subset-comm-props)]
    (mle/estimate-mixed-mult seed-objs-feats
                             other-props
                             (mx/mget subgraph-memb comm-idx))))

;; TODO: can we limit the number of vector ops by rounding some community prop values to zero?
;;       appears we would _not_ want to do this for object features, since those are tf-idf and possibly
;;       not as easy to determine if rounding to zero is not harmful to esimation.
(s/defn score-params :- Vec
  [groups-objs :- BinMat
   objs-feats :- Mat
   comms-props :- Mat]

  ;; (pmap (comp (partial reduce mx/add! (SparseIndexedVector/createLength (mx/column-count objs-feats)))
  ;;                                              ;;(fn [x] (println (mx/shape x)) x)
  ;;                                              (partial restricted-objs groups-objs objs-feats))
  ;;                                        (range (mx/row-count comms-props)))

  ;; (pmap (comp (partial reduce mx/add! (Vectorz/wrap (double-array (mx/column-count objs-feats) 0.0)))
  ;;                                              ;;(fn [x] (println (mx/shape x)) x)
  ;;                                              (partial restricted-objs groups-objs objs-feats))
  ;;                                        (range (mx/row-count comms-props)))
  (let [num-feats (mx/column-count objs-feats)]
    (sparse-indexed-vector (pmap log-multicat
                                 (mx/rows comms-props)
                                 ;; TODO: Make the map over rows in reduce mx/add! explicit
                                 (pmap (comp (fn [^SparseRowMatrix restrict-objs-feats]
                                               (try 
                                                 (reduce
                                                  (fn [sum x]
                                                    (mx/add! sum x))
                                                  (create-sparse-indexed-vector num-feats)
                                                  restrict-objs-feats)
                                                 (catch ArrayIndexOutOfBoundsException e
                                                   (let [out-file (->> (System/getProperty "user.dir")
                                                                       (java.io.File/createTempFile "oob-" ".clj"))]
                                                     (write-log (format "Logged data related to OOB exception at: %s." (.getAbsolutePath out-file)))
                                                     (spit out-file (for [^SparseIndexedVector row (.getRows restrict-objs-feats)
                                                                          :let [^Index index (get-private-field row "index")
                                                                                ^doubles data (get-private-field row "data")]]
                                                                      {:index (-> index .data seq)
                                                                       :data (-> data seq)})))  
                                                   (.printStackTrace e System/out)
                                                   (.flush System/out))))
                                             (partial restricted-objs groups-objs objs-feats))
                                       (range (mx/row-count comms-props)))))))


(s/defn alt-mixed-prop :- ProbVec
  "Create a mixture distribution given a set of objects, their membership weights,
  and community feature distributions."
  [groups-objs :- BinMat
   objs-groups :- BinMat   
   objs-feats :- Mat
   objs-membs :- Mat
   comms-props :- Mat
   comm-idx :- s/Int]
  (let [seed-objs (mx/get-row groups-objs comm-idx) ; objs part of the selected community seed subgraph
        subset-comms-props (selected-rows comms-props
                                          (apply b-or (selected-rows objs-groups seed-objs)))
        ;; TODO: use dense vectors here? add! is slow...
        seed-objs-feats (reduce mx/add!
                                (SparseIndexedVector/createLength (mx/column-count objs-feats))
                                (mx/rows (selected-rows objs-feats seed-objs)))
        
        ;; contribution proportion based on number of features per-object, N x 1
        obj-prop-num-obs (proportional (reduce mx/add!
                                               (SparseIndexedVector/createLength (mx/row-count objs-feats))
                                               (mx/columns (selected-rows objs-feats
                                                                          seed-objs))))
        ;; the complete phi vector
        subgraph-memb (->> (subgraph-membership (selected-rows objs-membs seed-objs)
                                                obj-prop-num-obs)
                           (round-to-zero! 1e-4)
                           proportional)]
    ;; These are SparseIndexedVector and SparseRowMatrix
    (mx/mmul subgraph-memb subset-comms-props)))

(s/defn m-step :- {:comms-props Mat
                    :changed s/Bool
                    :change-count s/Int
                    :changed-ids #{s/Int}}
  "Calculate new community topic distribution parameters.
  The groups-objs matrix can be used to specify which objects should
  represent a specific community. E.g., a clique vs. all objects
  with non-zero membership probabilty."
  [groups-objs :- BinMat
   objs-groups-possible :- BinMat
   objs-membs :- Mat
   objs-feats :- Mat
   comms-props :- Mat]
  ;; We want to score using all objects since we estimated comms with only seed group objects.
  ;; TODO: another conversion from rows to columns
  (let [scorer (partial score-params (sparse-row-matrix (mx/columns objs-groups-possible)) objs-feats)
        num-comms (mx/row-count comms-props)
        new-props (->> (sparse-row-matrix (pmap (partial estimate-comm-props
                                                         groups-objs
                                                         objs-groups-possible
                                                         objs-membs
                                                         objs-feats
                                                         comms-props)
                                                (range num-comms)))
                       (round-to-zero! 1e-4))
        old-mixed-props (sparse-row-matrix (pmap (partial alt-mixed-prop
                                                          groups-objs
                                                          objs-groups-possible
                                                          objs-feats
                                                          objs-membs
                                                          comms-props)
                                                 (range num-comms)))
        new-mixed-props (sparse-row-matrix (pmap (partial alt-mixed-prop
                                                          groups-objs
                                                          objs-groups-possible
                                                          objs-feats
                                                          objs-membs
                                                          new-props)
                                                     (range num-comms)))
        changed (atom false)
        change-count (atom 0)
        changed-ids (atom #{})]
    {:comms-props (sparse-row-matrix (pmap (fn [[old-score old-v] [new-score new-v] comm-idx]
                                             (if (>= old-score new-score)
                                               old-v
                                               (do (swap! changed not)
                                                   (swap! change-count inc)
                                                   (swap! changed-ids conj comm-idx)
                                                   new-v)))
                                           (map vector (scorer old-mixed-props) comms-props)
                                           (map vector (scorer new-mixed-props) new-props)
                                           (range num-comms)))
     :changed @changed
     :change-count @change-count
     :changed-ids @changed-ids}))

(s/defn run
  "Run the algorithm. Main calls this, so can you.
  Arguments:
  objs-feat-vals - observed feature values
  groups-objs - the seed objects of each group
  objs-groups - the subgraph groups an object may be affiliated with, does not change
  objs-membs - initial estimate of object membership to communities
  comms-props - initial estimate of community topic parameters
  max-iter - maximum number of iterations"
  [objs-feat-vals :- Mat
   groups-objs :- BinMat
   objs-groups :- BinMat
   objs-membs :- Mat
   comms-props :- Mat
   max-iter :- s/Int]
  (let [num-comms (mx/row-count comms-props)]
    (loop [membs objs-membs
           props comms-props
           iter-count 0
           last-changed-comms #{}
           all-changed-comms #{}]

      (write-log (format "Iteration: %d" iter-count))
      
      (if (= max-iter iter-count)
        (do (write-log "Maximum number of iterations reached.")
            {:objs-membs membs
             :comms-props props})
        (let [new-membs (e-step (partial restricted-comms objs-groups) objs-feat-vals props)
              {changed :changed
               change-count :change-count
               changed-ids :changed-ids
               new-props :comms-props} (m-step groups-objs objs-groups new-membs objs-feat-vals props)]
          (if (not changed)
            (do (write-log "Community distributions have converged.")
                {:objs-membs membs
                 :comms-props props})
            (do (write-log (format "Number of communities changed: %d" change-count))
                (write-log (format "Number of different communities from last iteration: %d" (count (set/difference changed-ids last-changed-comms))))
                (write-log (format "Number of first-time changed communities: %d" (count (set/difference changed-ids all-changed-comms))))
                (recur new-membs
                       new-props
                       (inc iter-count)
                       changed-ids
                       (set/union all-changed-comms changed-ids)))))))))

;; TODO: Add schemas for matrices to ensure dimensions are correct. For example, the objs-feat-vals matrix should have
;; the same number of columns as the comm-params matrix. I.e., we want a poor man's dependent types from predicate
;; schemas.
(defn -main
  [& args]
  (let [[opts args banner] (cli args
                                ["-h" "--help" "Print this help message and exit" :default false :flag true]
                                ["-o" "--out-dir" "Output directory" :default (System/getProperty "user.dir")]
                                ["-m" "--max-iter" "Max iterations" :default 10 :parse-fn #(Integer/parseInt %)])]
    (when (or (:help opts) (not= 1 (count args)))
      (println banner)
      (System/exit 0))

    (write-log "Starting up...")
    
    (let [input-path (nth args 0)
          max-iter (:max-iter opts)
          out-dir (:out-dir opts)
          objs-feat-vals (load-matrix (format "%s/objs-feat-vals.mm" input-path))
          ;; upper-bound, objs x comms
          candidate-comms (load-matrix (format "%s/upper-bound-objs-groups.mm" input-path))
          ;; lower-bound, comms x objs
          lower-bound-membs (load-matrix (format "%s/lower-bound-groups-objs.mm" input-path))
          _ (write-log "Loaded input.")
          initial-comms-props (->> (estimate-props objs-feat-vals lower-bound-membs)
                                   (round-to-zero! 1e-4))
          _ (write-log "Initial community proportions found.")
          ;; _ (write-log (format "Found %s non-zero communities." (->> initial-comms-props
          ;;                                                            mx/rows
          ;;                                                            (map mx/non-zero-count)
          ;;                                                            (filter (comp not zero?)) count)))
          ;; All the objects in at least one seed group
          seed-objs-idx (->> (mx/non-zero-indices lower-bound-membs)
                             (mapcat identity)
                             (distinct))
          num-objs (mx/column-count lower-bound-membs)
          initial-objs-membs (estimate-membs objs-feat-vals candidate-comms initial-comms-props (range num-objs))
          _ (write-log "Initial object memberships found.")
          {est-comms-props :comms-props
           est-objs-membs :objs-membs} (run
                                         objs-feat-vals
                                         lower-bound-membs
                                         candidate-comms
                                         initial-objs-membs
                                         initial-comms-props
                                         max-iter)
          _ (write-log "Estimating all objects membership...")
          num-objs (mx/row-count objs-feat-vals)
          all-est-objs-membs (estimate-membs objs-feat-vals candidate-comms est-comms-props (range num-objs))]
      
      (write-log "Saving output...")
      (save-matrix est-comms-props (format "%s/mle-comms-props.mm" out-dir))
      (save-matrix all-est-objs-membs (format "%s/mle-objs-membs.mm" out-dir)))))
