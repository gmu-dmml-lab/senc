(ns senc.util-test
  (:require [clojure.test :refer :all]
            [senc.util :refer :all]
            [clojure.core.matrix :as mx]
            [schema.test]))

(use-fixtures :once schema.test/validate-schemas)

(mx/set-current-implementation :vectorz)

(defn close? [^double delta-thresh ^double expected ^double actual]
  (< (Math/abs (- expected actual)) delta-thresh))

(defn close-seq? [^double delta-thresh expected actual]
  (every? (fn [[^double e ^double a]] (< (Math/abs (- e a)) delta-thresh))
          (map vector expected actual)))
