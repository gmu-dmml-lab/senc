(ns senc.util
  (:require [clojure.core.matrix :as mx]
            [schema.core :as s]
            [munge.matrix :as mmx])
  (:use [munge.schema :only [Vec Mat ProbVec BinVec]])
  (:import [mikera.matrixx.impl SparseRowMatrix]))

(set! *warn-on-reflection* true)

(s/defn safe-log :- s/Num
  "A log function that returns a very small number for input value zero.
  Simplifies working with probability vectors where zero-value
  elements can be ignored."
  [x :- s/Num]
  (if (zero? x)
    1e-10
    (Math/log x)))

(s/defn log-multicat :- s/Num
  "This is a bastardized not-quite-legit log PMF. It's essentially
  a multinomial without the multinomial coefficient. This is useful
  for computing relative change in likelihood for a distribution
  we're trying to estimate."
  [p :- Vec
   xs :- Vec]
  ;; TODO: using 0's instead of 1e-10 would speed up inner product
  (mx/mget ^mikera.vectorz.AScalar (mx/inner-product xs (mx/emap safe-log p))))
