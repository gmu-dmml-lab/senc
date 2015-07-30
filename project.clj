(defproject senc "0.1.1"
  :description "Seeded estimation of network communities. Implements the method proposed in \"Finding Community Topics and Membership in Graphs\"."
  :url "https://github.com/gmu-dmml-lab/senc"
  :license {:name "MIT"
            :url "http://opensource.org/licenses/MIT"}
  :dependencies [[org.clojure/clojure "1.7.0"]
                 [org.clojure/tools.cli "0.3.1"]
                 [prismatic/schema "0.4.3"]
                 [notion "0.1.1"]
                 [munge "0.1.1"]
                 ;;[net.mikera/core.matrix "0.33.3-SNAPSHOT"]
                 ;;[net.mikera/vectorz-clj "0.29.1-SNAPSHOT"]
                 [net.mikera/core.matrix "0.36.1"]
                 [net.mikera/vectorz-clj "0.30.1"]]
  :main senc.core
  :profiles {:common {:jvm-opts ["-Xmx14g" "-server" "-XX:+UseConcMarkSweepGC" "-XX:+TieredCompilation"]}
             :prod [:common]}
  :global-vars {*warn-on-reflection* true}
  :test-selectors {:default (complement :slow)
                   :slow :slow
                   :all (constantly true)})
