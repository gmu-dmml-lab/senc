(ns senc.estimate.mle
  (:require [clojure.core.matrix :as mx]
            [schema.core :as s])
  (:use [munge.schema :only [Vec Mat ProbVec]]
        [munge.matrix :only [proportional select-rows]]))

(s/defn solve-linear :- s/Num
  "Solve a simple linear equation of the form y = ax + b."
  [y :- s/Num
   a :- s/Num
   b :- s/Num]
  (/ (- y b) a))

(s/defn estimate-mixed-mult :- ProbVec
  "We are solving a series of linear equations, 
  one for each categorical dimension.

  The observations, other community proportions, and 
  mixing weight of the community being estimated are 
  provided."
  [obs :- Vec
   ;; TODO: when other-props isn't ProbVec is it all zeros? make a ZeroVec?
   other-props :- (s/either ProbVec Vec)
   mix-weight :- s/Num]
  ;; final proportional needed when there were negative values from subtracting the
  ;; proportion of observed feature values by the mixed other community parameters
  (->> (-> (proportional obs)
           (mx/sub! other-props)
           (mx/div! mix-weight))
       (mx/emap! (partial max 0))
       (proportional)))

(s/defn estimate-comm :- Vec
  [select-objs :- (s/either (s/=> [s/Int] s/Int) [[s/Int]])
   obj-feats :- Mat
   obj-memb :- Mat
   comm-props :- Mat
   comm-idx :- s/Int]
  (let [obj-idxs (select-objs comm-idx)
        comb-memb (proportional (reduce mx/add (select-rows obj-memb obj-idxs)))
        comb-feats (reduce mx/add (select-rows obj-feats obj-idxs))
        other-props (mx/mmul (mx/mset comb-memb comm-idx 0) comm-props)
        mix-weight (mx/mget comb-memb comm-idx)]
    (estimate-mixed-mult comb-feats other-props mix-weight)))

(s/defn estimate-mix-from-terms :- Vec
  "Estimate the mixture of communities for an object."
  [obj-feats :- Mat
   obj-memb :- Mat
   comm-props :- Mat
   obj-idx :- s/Int]
  (assert false "Unimplemented."))
