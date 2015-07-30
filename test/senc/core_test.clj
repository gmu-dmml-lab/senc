(ns senc.core-test
  (:require [clojure.test :refer :all]
            [senc.util-test :refer [close? close-seq?]]
            [senc.core :refer :all]
            [munge.matrix :refer [selected-rows sparse-indexed-vector sparse-matrix]]
            [clojure.core.matrix :as mx]
            [schema.test]
            [clojure.core.matrix.impl.pprint :as mpp :refer [pm]]))

(use-fixtures :once schema.test/validate-schemas)

(mx/set-current-implementation :vectorz)

(deftest estimate-memb-test
  (testing "object is member of only one group"
    (is (close-seq? 1e-3
                    [1.0 0.0]
                    (estimate-memb (sparse-indexed-vector [10 0 0])
                                   (selected-rows (sparse-matrix 2 3 [[0.4 0.4 0.2]
                                                                      [0.1 0.8 0.1]])
                                                  (sparse-indexed-vector [1 0]))))))
  (testing "basic membership"
    (is (close-seq? 1e-4
                    [0.3333 0.6666]
                    (estimate-memb (sparse-indexed-vector [10 0 20 0 30])
                                   (sparse-matrix 2 4[[0.6 0.3 0.1 0 0]
                                                      [0 0 0.1 0.3 0.6]]))))
    (is (close-seq? 1e-4
                    [0 0.477 0.523 0 0]
                    (estimate-memb (sparse-indexed-vector [51 0 20 8 21])
                                   (selected-rows (sparse-matrix 5 5 [[0.7 0.2 0.1 0 0]
                                                                      [0.7 0 0.3 0 0]
                                                                      [0.3 0 0.2 0.2 0.3]
                                                                      [0.4 0.5 0.1 0 0]
                                                                      [0.1 0.1 0.2 0.4 0.2]])
                                                  (sparse-indexed-vector [0 1 1 0 0]))))))

  (testing "disjoint features"
    (is (mx/equals [0.5 0.5]
                   (estimate-memb (sparse-indexed-vector [10 5 10 5])
                                  (sparse-matrix 2 3 [[0.9 0.1 0 0]
                                                      [0 0 0.9 0.1]]))))
    (testing "feature count changes the mix"
      (is (mx/equals [0.75 0.25]
                     (estimate-memb (sparse-indexed-vector [10 5 3 2])
                                    (sparse-matrix 2 3 [[0.9 0.1 0 0]
                                                        [0 0 0.9 0.1]]))))))
 
  (testing "merging membership vectors"
    ;; TODO: test if estimating a group of objects directly results
    ;; in same final membership weight as doing individually and merging
    ;; as we do in m-step
    ))

(deftest estimate-membs-test
  (testing "basic"
    (is (mx/equals (mx/matrix [[1.0 0.0]
                               [0.0 0.0]
                               [0.0 1.0]])
                   (estimate-membs (mx/matrix [[30 20 10]
                                               [10 10 10]
                                               [20 50 10]])
                                   (mx/matrix [[1 0]
                                               [1 1]
                                               [0 1]])
                                   (mx/matrix [[0.5 0.3 0.2]
                                               [0.25 0.6 0.15]])
                                   [0 2])
                   1e-2))))

(deftest restricted-comms-test
  (is (mx/equals (mx/matrix [[0.5 0.3 0.2]
                             [0.7 0.1 0.2]
                             [0 0 0]])
                 (restricted-comms (mx/matrix [[1 0 0]
                                               [1 1 0]
                                               [1 0 0]])
                                   (mx/matrix [[0.5 0.3 0.2]
                                               [0.7 0.1 0.2]
                                               [0.4 0.4 0.2]])
                                   #{0 1})))
  (is (mx/equals (mx/matrix [[0.5 0.2 0.3]
                             [0.0 0.0 0.0]])
                 (restricted-comms (mx/matrix [[1 0]])
                                   (mx/matrix [[0.5 0.2 0.3]
                                               [0.1 0.1 0.8]])
                                   #{0}))))

(deftest e-step-test
  (testing "trivial case of single membership"
    (is (mx/equals (mx/matrix [[1 0 0]
                               [0 1 0]
                               [0 0 1]])
                   (e-step (partial restricted-comms (mx/matrix [[1 0 0]
                                                                 [0 1 0]
                                                                 [0 0 1]]))
                           (mx/matrix [[5 2 0 0]
                                       [1 3 0 0]
                                       [0 0 3 10]])
                           (mx/matrix [[0.4 0.2 0.2 0.2]
                                       [0.2 0.4 0.2 0.0]
                                       [0.0 0.0 0.3 0.7]])))))
  (testing "mixture membership"
    (is (close-seq? 1e-2
                    (concat [0.612 0.388 0.0]
                            [0.433 0.566 0.0]
                            [0.521 0.0 0.479]
                            [0.0 0.41 0.59])
                    (mx/eseq (e-step (partial restricted-comms (mx/matrix [[1 1 0]
                                                                           [1 1 0]
                                                                           [1 0 1]
                                                                           [0 1 1]]))
                                     (mx/matrix [[5 2 1 1]
                                                 [1 3 1 0]
                                                 [5 2 3 10]
                                                 [5 2 3 10]])
                                     (mx/matrix [[0.4 0.2 0.2 0.2]
                                                 [0.2 0.4 0.2 0.0]
                                                 [0.0 0.0 0.3 0.7]])))))))

(deftest subgraph-membership-test
  (testing "trivial case"
    (is (close-seq? 1e-2
                    [0.84 0.12 0.04]
                    (subgraph-membership (mx/matrix [[1.0 0.0 0.0]
                                                     [0.6 0.3 0.1]])
                                         (mx/matrix [0.6
                                                     0.4]))))))

(deftest mixed-props-test
  (testing "trivial case when objects only belong to community being estimated"
    ;; For comm-idx 0
    (is (close-seq? 1e-2
                    [0.0 0.0 0.0]
                    (mixed-props (mx/matrix [0.0 0.0 0.0])
                                 (mx/matrix [[0.6 0.4 0.0]
                                             [0.3 0.7 0.0]
                                             [0.9 0.0 0.1]])))))
  (testing "simple case"
    ;; For comm-idx 0
    (is (close-seq? 1e-2
                    [0.30 0.28 0.02]
                    (mixed-props (mx/matrix [0.0 0.4 0.2])
                                 (mx/matrix [[0.6 0.4 0.0]
                                             [0.3 0.7 0.0]
                                             [0.9 0.0 0.1]]))))
    ;; For comm-idx 1
    (is (close-seq? 1e-2
                    [0.42 0.16 0.02]
                    (mixed-props (mx/matrix [0.4 0.0 0.2])
                                 (mx/matrix [[0.6 0.4 0.0]
                                             [0.3 0.7 0.0]
                                             [0.9 0.0 0.1]]))))))

(deftest estimate-comm-props-test
  (testing "trivial case of estimating community parameters"
    (is (close-seq? 1e-2
                    [0.66 0.33 0 0]
                    (mx/eseq (estimate-comm-props (mx/matrix [[1 0]
                                                              [0 1]])
                                                  (mx/matrix [[1 0]
                                                              [0 1]])
                                                  (mx/matrix [[1.0 0]
                                                              [0 1.0]])
                                                  (mx/matrix [[10 5 0 0]
                                                              [0 0 10 10]])
                                                  (mx/matrix [[0.7 0.3 0 0]
                                                              [0 0 0.6 0.4]])
                                                  0)))))
  (testing "case of estimating community parameters where result sums greater 
            than one b/c some features occurred less than expected and 
            not were overaccoutned for by other community distributions"
    (is (close-seq? 1e-2
                    [0.69 0.31 0.0 0.0]
                    (mx/eseq (estimate-comm-props (mx/matrix [[1 0 1]
                                                              [0 0 0]
                                                              [0 1 1]])
                                                  (mx/matrix [[1 0 0]
                                                              [0 0 1]
                                                              [1 0 1]])
                                                  (mx/matrix [[1.0 0 0]
                                                              [0 1.0 0]
                                                              [0.3 0 0.7]])
                                                  (mx/matrix [[10 5 0 0]
                                                              [0 0 10 10]
                                                              [6 2 4 10]])
                                                  (mx/matrix [[0.7 0.3 0 0]
                                                              [0 0 0.6 0.4]
                                                              [0 0 0.3 0.7]])
                                                  0))))))

(deftest alt-mixed-prop-test
  (testing "basic case"
    (is (close-seq? 1e-2
                    [0.5 0.3 0.2 0.0]
                    (alt-mixed-prop (mx/matrix [[1 0 0]
                                                [0 1 0]
                                                [0 0 1]])
                                    (mx/matrix [[1 0 0]
                                                [0 1 0]
                                                [0 0 1]])
                                    (mx/matrix [[50 25 25 0]
                                                [10 10 10 10]
                                                [10 20 30 40]])
                                    (mx/matrix [[1.0 0.0 0.0]
                                                [0.0 1.0 0.0]
                                                [0.0 0.0 1.0]])
                                    (mx/matrix [[0.5 0.3 0.2 0.0]
                                                [0.25 0.25 0.25 0.25]
                                                [0.1 0.2 0.3 0.4]])
                                    0))))
  (testing "simple mixture"
    (is (close-seq? 1e-2
                    [0.375 0.275 0.225 0.125]
                    (alt-mixed-prop (mx/matrix [[1 0 0]
                                                [1 1 0]
                                                [0 0 1]])
                                    (mx/matrix [[1 1 0]
                                                [0 1 0]
                                                [0 0 1]])
                                    (mx/matrix [[50 25 25 0]
                                                [10 10 10 10]
                                                [10 20 30 40]])
                                    (mx/matrix [[0.5 0.5 0.0]
                                                [0.0 1.0 0.0]
                                                [0.0 0.0 1.0]])
                                    (mx/matrix [[0.5 0.3 0.2 0.0]
                                                [0.25 0.25 0.25 0.25]
                                                [0.1 0.2 0.3 0.4]])
                                    0)))))

;; num object groups = 5
;; num comms = 5
;; num features = 6
(deftest ^:slow run-uniform-test
  (testing "simple set of communities and objects"
    (let [feat-vals (mx/matrix [[73  8 19  0  0]
                                [72  1 27  0  0]
                                [69  1 30  0  0]
                                [65  9 26  0  0]
                                [68 22 10  0  0]
                                [72 14 14  0  0]
                                [68  0 32  0  0]
                                [74 16 10  0  0]
                                [66  5 29  0  0]
                                [64 22 14  0  0]
                                [71 11 18  0  0]
                                [75  5 20  0  0]
                                [66 17 17  0  0]
                                [76 13 11  0  0]
                                [66  1 33  0  0]
                                [77  7 16  0  0]
                                [70 11 19  0  0]
                                [72 14 14  0  0]
                                [69 14 17  0  0]
                                [66 12 22  0  0]
                                [68 13 19  0  0]
                                [68  0 32  0  0]
                                [67 10 23  0  0]
                                [70  6 24  0  0]
                                [76  1 23  0  0]
                                [67 19 14  0  0]
                                [64 16 20  0  0]
                                [62 22 16  0  0]
                                [74 12 14  0  0]
                                [65  0 35  0  0]
                                [65  5 30  0  0]
                                [65 15 20  0  0]
                                [71 16 13  0  0]
                                [68 15 17  0  0]
                                [76  1 23  0  0]
                                [68  9 23  0  0]
                                [65 21 14  0  0]
                                [70  6 24  0  0]
                                [72  8 20  0  0]
                                [70  0 30  0  0]
                                [80  6 14  0  0]
                                [77 13 10  0  0]
                                [74  6 20  0  0]
                                [73  9 18  0  0]
                                [63  5 32  0  0]
                                [76  0 24  0  0]
                                [67 16 17  0  0]
                                [66 21 13  0  0]
                                [76  0 24  0  0]
                                [67 14 19  0  0]
                                [71  8 21  0  0]
                                [71 18 11  0  0]
                                [71 11 18  0  0]
                                [70 22  8  0  0]
                                [78 11 11  0  0]
                                [75  6 19  0  0]
                                [71 13 16  0  0]
                                [72 16 12  0  0]
                                [73 15 12  0  0]
                                [75  6 19  0  0]
                                [32  3 20 17 28]
                                [66  4 25  2  3]
                                [56 11 17  6 10]
                                [56 11 15  8 10]
                                [43  3 21 11 22]
                                [27  4 17 23 29]
                                [55  4 21  9 11]
                                [50  4 24  6 16]
                                [55 17 25  0  3]
                                [45  5 20 13 17]
                                [75  9 16  0  0]
                                [51  7 17 11 14]
                                [37  2 31 14 16]
                                [64 10 18  3  5]
                                [64  0 22  3 11]
                                [55  1 27  7 10]
                                [52  6 15 14 13]
                                [61  7 21  5  6]
                                [53  5 22  4 16]
                                [37  4 20 13 26]
                                [63 14 11  8  4]
                                [75 13 11  1  0]
                                [63  1 29  3  4]
                                [69 14 16  0  1]
                                [34  0 24 15 27]
                                [47  8 17 11 17]
                                [63 14 12  4  7]
                                [68  5 26  1  0]
                                [49 11 16  8 16]
                                [68  4 22  2  4]
                                [37  0 21 10 32]
                                [61  0 28  7  4]
                                [60  0 30  2  8]
                                [57  0 30  3 10]
                                [47  0 28 13 12]
                                [70  0 29  0  1]
                                [63  0 36  0  1]
                                [49  0 21 12 18]
                                [51  0 26  4 19]
                                [71  0 29  0  0]
                                [51 32 14  2  1]
                                [48 24 18  7  3]
                                [45 24 16  7  8]
                                [34 57  5  0  4]
                                [43 25 11 15  6]
                                [61  8 20  8  3]
                                [27 10 21 28 14]
                                [54 16 25  3  2]])
          obj-groups (mx/matrix (concat (repeat 60 [1 1 0 0 0])
                                        (repeat 30 [1 1 1 0 0])
                                        (repeat 10 [0 1 1 0 0])
                                        (repeat 5 [1 0 0 1 1])
                                        (repeat 3 [0 1 0 1 1])))
          ;; NOTE: not usually the case this is the transpose
          groups-objs (mx/transpose obj-groups)
          cores (mx/matrix [(concat (repeat 60 1) (repeat 48 0))
                            (concat (repeat 60 0) (repeat 30 1) (repeat 18 0))
                            (concat (repeat 90 0) (repeat 10 1) (repeat 8 0))
                            (concat (repeat 100 0) (repeat 5 1) (repeat 3 0))
                            (concat (repeat 105 0) (repeat 3 1))])
          init-comm-props (estimate-props feat-vals cores)
          init-obj-membs (mx/matrix (map (fn [obj-feat-vals obj-idx]
                                           (estimate-memb obj-feat-vals
                                                          (restricted-comms obj-groups
                                                                            init-comm-props
                                                                            (hash-set obj-idx))))
                                         feat-vals
                                         (range (mx/row-count feat-vals))))]
      (is (close-seq? 1e-2
                      '(0.703 0.138 0.159 0.0 0.0 0.605 0.086 0.190 0.049 0.069 0.462 0.001 0.253 0.103 0.178 0.272 0.505 0.079 0.099 0.044 0.289 0.309 0.167 0.170 0.063)
                      (-> (run
                            feat-vals
                            groups-objs
                            obj-groups
                            init-obj-membs
                            init-comm-props
                            20)
                          :comms-props
                          mx/eseq))))))  

;; num object groups = 5
;; num comms = 5
;; num features = 6
;; generated with seed at 1001
(deftest ^:slow run-new-uniform-test
  (testing "simple set of communities and objects"
    (let [feat-vals (mx/matrix [[73 19 8 0 0]   
                                [66 17 17 0 0]   
                                [80 4 16 0 0]   
                                [72 7 21 0 0]   
                                [70 10 20 0 0]   
                                [79 5 16 0 0]   
                                [63 18 19 0 0]   
                                [68 3 29 0 0]   
                                [70 5 25 0 0]   
                                [65 15 20 0 0]   
                                [71 17 12 0 0]   
                                [72 5 23 0 0]   
                                [77 2 21 0 0]   
                                [76 12 12 0 0]   
                                [70 15 15 0 0]   
                                [67 11 22 0 0]   
                                [70 0 30 0 0]   
                                [81 9 10 0 0]   
                                [71 18 11 0 0]   
                                [71 8 21 0 0]   
                                [65 13 22 0 0]   
                                [69 5 26 0 0]   
                                [66 21 13 0 0]   
                                [66 7 27 0 0]   
                                [69 18 13 0 0]   
                                [75 15 10 0 0]   
                                [69 5 26 0 0]   
                                [67 27 6 0 0]   
                                [79 2 19 0 0]   
                                [68 18 14 0 0]   
                                [74 5 21 0 0]   
                                [69 11 20 0 0]   
                                [63 18 19 0 0]   
                                [76 3 21 0 0]   
                                [70 13 17 0 0]   
                                [69 5 26 0 0]   
                                [73 1 26 0 0]   
                                [70 7 23 0 0]   
                                [74 8 18 0 0]   
                                [77 8 15 0 0]   
                                [63 4 33 0 0]   
                                [69 7 24 0 0]   
                                [73 21 6 0 0]   
                                [64 11 25 0 0]   
                                [67 12 21 0 0]   
                                [70 14 16 0 0]   
                                [75 9 16 0 0]   
                                [77 10 13 0 0]   
                                [70 0 30 0 0]   
                                [72 6 22 0 0]   
                                [73 8 19 0 0]   
                                [76 2 22 0 0]   
                                [74 19 7 0 0]   
                                [63 22 15 0 0]   
                                [61 27 12 0 0]   
                                [64 22 14 0 0]   
                                [66 8 26 0 0]   
                                [73 17 10 0 0]   
                                [73 6 21 0 0]   
                                [70 6 24 0 0]   
                                [66 0 29 1 4]   
                                [53 6 20 12 9]   
                                [50 3 18 8 21]  
                                [48 5 29 9 9]   
                                [74 7 15 3 1]   
                                [59 11 21 6 3]   
                                [68 4 14 7 7]   
                                [70 0 23 4 3]   
                                [48 2 25 9 16]  
                                [59 9 13 4 15]  
                                [41 0 26 15 18]  
                                [62 7 27 2 2]   
                                [43 6 17 14 20]  
                                [46 3 19 15 17]  
                                [55 5 20 8 12]  
                                [37 1 27 12 23]  
                                [64 2 21 6 7]   
                                [63 18 16 0 3]   
                                [50 9 14 9 18]  
                                [62 15 14 2 7]   
                                [47 1 22 12 18]  
                                [69 3 25 2 1]   
                                [58 13 14 8 7]   
                                [63 7 20 4 6]   
                                [64 1 22 9 4]   
                                [69 11 17 0 3]   
                                [49 9 26 11 5]   
                                [63 12 21 1 3]   
                                [72 0 25 1 2]   
                                [44 2 16 17 21]  
                                [51 0 20 8 21]  
                                [58 0 28 7 7]   
                                [34 0 32 10 24]  
                                [57 0 24 6 13]  
                                [47 0 14 13 26]  
                                [69 0 28 0 3]   
                                [43 0 19 15 23]  
                                [42 0 26 10 22]  
                                [66 0 29 2 3]   
                                [64 0 30 2 4]   
                                [65 28 7 0 0]   
                                [47 31 12 5 5]   
                                [29 28 11 18 14]  
                                [26 26 18 21 9]   
                                [30 14 17 27 12]  
                                [20 16 17 28 19]  
                                [34 25 15 16 10]  
                                [34 11 21 23 11]])
          obj-groups (mx/matrix (concat (repeat 60 [1 1 0 0 0])
                                        (repeat 30 [1 1 1 0 0])
                                        (repeat 10 [0 1 1 0 0])
                                        (repeat 5 [1 0 0 1 1])
                                        (repeat 3 [0 1 0 1 1])))
          cores (mx/matrix [(concat (repeat 60 1) (repeat 48 0))
                            (concat (repeat 60 0) (repeat 30 1) (repeat 18 0))
                            (concat (repeat 90 0) (repeat 10 1) (repeat 8 0))
                            (concat (repeat 100 0) (repeat 5 1) (repeat 3 0))
                            (concat (repeat 105 0) (repeat 3 1))])
          ;; NOTE: not usually the case this is the transpose
          groups-objs (mx/transpose obj-groups)
          init-comm-props (estimate-props feat-vals cores)
          init-obj-membs (mx/matrix (map (fn [obj-feat-vals obj-idx]
                                           (estimate-memb obj-feat-vals
                                                          (restricted-comms obj-groups
                                                                            init-comm-props
                                                                            (hash-set obj-idx))))
                                         feat-vals
                                         (range (mx/row-count feat-vals))))
          result (run
                   feat-vals
                   groups-objs
                   obj-groups
                   init-obj-membs
                   init-comm-props
                   50)]
      (println (pm (:obj-membs result)))
      (is (close-seq? 1e-2
                      '(0.706 0.107 0.188 0.0 0.0 0.572 0.057 0.205 0.070 0.095 0.531 0.0 0.25 0.073 0.146 0.394 0.254 0.13 0.142 0.08 0.293 0.173 0.176 0.223 0.133)
                      (-> result
                          :comms-props
                          mx/eseq))))))


(deftest ^:slow run-select-heavy-test
  (testing "simple set of communities and objects"
    (let [feat-vals (mx/matrix [[76 15 9  0  0]   
                                [66 23 11 0  0]   
                                [70 21 9  0  0]   
                                [77 12 11 0  0]   
                                [73 13 14 0  0]   
                                [70 15 15 0  0]   
                                [77 11 12 0  0]   
                                [62 19 19 0  0]   
                                [71 19 10 0  0]   
                                [70 21 9  0  0]   
                                [66 20 14 0  0]   
                                [74 17 9  0  0]   
                                [70 23 7  0  0]   
                                [79 10 11 0  0]   
                                [71 14 15 0  0]   
                                [71 19 10 0  0]   
                                [70 20 10 0  0]   
                                [66 24 10 0  0]   
                                [81 9  10 0  0]   
                                [70 17 13 0  0]   
                                [70 18 12 0  0]   
                                [68 19 13 0  0]   
                                [71 17 12 0  0]   
                                [61 23 16 0  0]   
                                [68 24 8  0  0]   
                                [73 8  19 0  0]   
                                [70 19 11 0  0]   
                                [72 17 11 0  0]   
                                [68 25 7  0  0]   
                                [75 21 4  0  0]   
                                [68 23 9  0  0]   
                                [74 17 9  0  0]   
                                [70 21 9  0  0]   
                                [67 15 18 0  0]   
                                [73 18 9  0  0]   
                                [71 20 9  0  0]   
                                [71 17 12 0  0]   
                                [71 20 9  0  0]   
                                [70 17 13 0  0]   
                                [74 8  18 0  0]   
                                [74 14 12 0  0]   
                                [66 24 10 0  0]   
                                [71 19 10 0  0]   
                                [68 23 9  0  0]   
                                [62 23 15 0  0]   
                                [69 18 13 0  0]   
                                [71 16 13 0  0]   
                                [74 14 12 0  0]   
                                [79 13 8  0  0]   
                                [70 25 5  0  0]   
                                [71 19 10 0  0]   
                                [75 14 11 0  0]   
                                [75 15 10 0  0]   
                                [72 16 12 0  0]   
                                [61 24 15 0  0]   
                                [67 28 5  0  0]   
                                [59 23 18 0  0]   
                                [67 19 14 0  0]   
                                [73 16 11 0  0]   
                                [75 14 11 0  0]   
                                [60 4  24 3  9]   
                                [65 1  32 1  1]   
                                [66 3  30 0  1]   
                                [71 0  28 1  0]   
                                [60 1  38 0  1]   
                                [72 0  21 4  3]   
                                [76 0  21 1  2]   
                                [63 1  31 4  1]   
                                [70 0  26 3  1]   
                                [67 0  26 5  2]   
                                [70 1  26 1  2]   
                                [74 0  25 0  1]   
                                [54 4  28 5  9]   
                                [60 2  33 3  2]   
                                [70 1  28 0  1]   
                                [68 2  27 1  2]   
                                [66 5  26 1  2]   
                                [77 1  21 1  0]   
                                [63 1  27 4  5]   
                                [66 0  31 1  2]   
                                [71 0  21 5  3]   
                                [61 0  35 1  3]   
                                [64 1  30 2  3]   
                                [69 1  29 1  0]   
                                [74 8  17 0  1]   
                                [64 1  33 1  1]   
                                [63 4  29 2  2]   
                                [69 0  30 1  0]   
                                [63 1  36 0  0]   
                                [70 6  23 1  0]   
                                [33 0  21 15 31]  
                                [35 0  21 21 23]  
                                [21 0  27 21 31]  
                                [40 0  22 7  31]  
                                [33 0  30 15 22]  
                                [27 0  25 20 28]  
                                [28 0  29 13 30]  
                                [28 0  26 18 28]  
                                [39 0  21 10 30]  
                                [39 0  18 18 25]  
                                [36 42 16 4  2]   
                                [40 54 5  0  1]   
                                [35 44 15 3  3]   
                                [32 56 9  2  1]   
                                [35 40 16 7  2]   
                                [22 16 25 27 10]  
                                [23 15 15 35 12]  
                                [13 19 21 26 21]])
          obj-groups (mx/matrix (concat (repeat 60 [1 1 0 0 0])
                                        (repeat 30 [1 1 1 0 0])
                                        (repeat 10 [0 1 1 0 0])
                                        (repeat 5 [1 0 0 1 1])
                                        (repeat 3 [0 1 0 1 1])))
          ;; NOTE: not usually the case this is the transpose
          groups-objs (mx/transpose obj-groups)
          cores (mx/matrix [(concat (repeat 60 1) (repeat 48 0))
                            (concat (repeat 60 0) (repeat 30 1) (repeat 18 0))
                            (concat (repeat 90 0) (repeat 10 1) (repeat 8 0))
                            (concat (repeat 100 0) (repeat 5 1) (repeat 3 0))
                            (concat (repeat 105 0) (repeat 3 1))])
          init-comm-props (estimate-props feat-vals cores)
          init-obj-membs (mx/matrix (map (fn [obj-idx]
                                           (estimate-memb (mx/get-row feat-vals obj-idx)
                                                          (restricted-comms obj-groups
                                                                            init-comm-props
                                                                            (hash-set obj-idx))))
                                         (range (mx/row-count feat-vals))))]
      (is (close-seq? 1e-2
                      '(0.706 0.181 0.113 0.0 0.0 0.669 0.016 0.277 0.018 0.02 0.323 0.0 0.24 0.158 0.279 0.356 0.472 0.122 0.032 0.018 0.193 0.166 0.203 0.293 0.143)
                      (-> (run
                            feat-vals
                            groups-objs
                            obj-groups
                            init-obj-membs
                            init-comm-props
                            50)
                          :comms-props
                          mx/eseq))))))
