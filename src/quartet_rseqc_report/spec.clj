(ns quartet-rseqc-report.spec
  (:require [clojure.spec.alpha :as s]
            [spec-tools.core :as st]))

(s/def ::name
  (st/spec
   {:spec                string?
    :type                :string
    :description         "The name of the report"
    :swagger/default     ""
    :reason              "Not a valid report name"}))

(s/def ::description
  (st/spec
   {:spec                string?
    :type                :string
    :description         "Description of the report"
    :swagger/default     ""
    :reason              "Not a valid description."}))

(s/def ::filepath
  (st/spec
   {:spec                (s/and string? #(re-matches #"^oss:\/\/[a-zA-Z0-9_\/]+.*" %))
    :type                :string
    :description         "File path for covertor."
    :swagger/default     nil
    :reason              "The filepath must be string."}))

(s/def ::group
  (st/spec
   {:spec                string?
    :type                :string
    :description         "A group name which is matched with library."
    :swagger/default     []
    :reason              "The group must a string."}))

(s/def ::library
  (st/spec
   {:spec                string?
    :type                :string
    :description         "A library name."
    :swagger/default     []
    :reason              "The library must a string."}))

(s/def ::sample
  (st/spec
   {:spec                string?
    :type                :string
    :description         "A sample name."
    :swagger/default     []
    :reason              "The sample name must a string."}))

(s/def ::metadat-item
  (s/keys :req-un [::library
                   ::group
                   ::sample]))

(s/def ::metadata
  (s/coll-of ::metadat-item))

(def quartet-rseqc-report-params-body
  "A spec for the body parameters."
  (s/keys :req-un [::name ::filepath ::metadata]
          :opt-un [::description]))
