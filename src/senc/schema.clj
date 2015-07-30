(ns senc.schema
  (:require [clojure.core.matrix :as mx]
            [schema.core :as s]))

(def Vec
  "Any vector"
  (s/pred mx/vec? 'vec?))

(def Mat
  "Any matrix"
  (s/pred mx/matrix? 'matrix?))

(defn prob-vec?
  [x]
  (and (mx/vec? x)
       (> 0.001
          (Math/abs (double (- 1 (mx/esum x)))))))
(def ProbVec
  "Probabilty vector"
  (s/pred prob-vec? 'prob-vec?))

(defn bin-vec?
  [x]
  (and (mx/vec? x)
       (every? #(or (== 1 %) (== 0 %)) (mx/eseq x))))
(def BinVec
  "Binary vector"
  (s/pred bin-vec? 'bin-vec?))

(defn bin-mat?
  [x]
  (and (mx/matrix? x)
       (every? #(or (== 1 %) (== 0 %)) (mx/eseq x))))
(def BinMat
  "Binary matrix"
  (s/pred bin-mat? 'bin-mat?))

(s/defrecord Obj
    [id :- Integer
     feat-counts :- [Double]])

(s/defn obj
  [id :- Integer
   feat-counts :- [Double]]
  (Obj. id feat-counts))

(s/defrecord Community
    [id :- Integer
     feat-props :- [Double]
     members :- [senc.schema.Obj]])

(s/defn community :- Community
  [id :-  Integer
   feat-props :- [Double]
   members :- [senc.schema.Obj]]
  (Community. id feat-props members))
