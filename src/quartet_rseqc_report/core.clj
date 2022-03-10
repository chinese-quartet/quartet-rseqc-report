(ns quartet-rseqc-report.core
  (:require [tservice-core.tasks.http :as http-task]
            [quartet-rseqc-report.spec :as spec]
            [quartet-rseqc-report.task :as task]))

(def metadata
  (http-task/make-routes "quartet-rseqc-report" :ReportPlugin
                         {:method-type :post
                          :endpoint "quartet-rseqc-report"
                          :summary "Generate the QC Report for Quartet Proteomics data."
                          :body-schema spec/quartet-rseqc-report-params-body
                          :response-schema any?
                          :handler task/post-handler}))

(def events-init task/events-init)
