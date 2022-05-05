(ns quartet-rseqc-report.task
  (:require [quartet-rseqc-report.rseqc :as rseqc]
            [local-fs.core :as fs-lib]
            [tservice-core.plugins.env :refer [get-workdir make-remote-link add-env-to-path create-task! update-task!]]
            [tservice-core.plugins.util :as util]
            [clojure.data.json :as json]
            [clojure.string :as clj-str]
            [clojure.tools.logging :as log]
            [tservice-core.tasks.async :refer [publish-event! make-events-init]]
            [quartet-rseqc-report.version :as v]))

(defn date
  []
  (.format (java.text.SimpleDateFormat. "yyyy-MM-dd")
           (new java.util.Date)))

(defn update-process!
  [^String task-id ^Integer percentage]
  (let [record (cond
                 (= percentage 100) {:status "Finished"
                                     :percentage 100
                                     :finished_time (util/time->int (util/now))}
                 (= percentage -1) {:status "Failed"
                                    :finished_time (util/time->int (util/now))}
                 :else {:percentage percentage})
        record (merge {:id task-id} record)]
    (update-task! record)))

(defn update-log-process!
  "Update message into log file and process into database."
  [log-path coll task-id process]
  (spit log-path (json/write-str coll))
  (update-process! task-id process))

(defn post-handler
  [{:keys [body owner plugin-context uuid workdir]
    :as payload}]
  (log/info (format "Create a report with %s" payload))
  (let [{:keys [name filepath description metadata]
         :or {description (format "Quality control report for %s" name)}} body
        payload (merge {:description description} (:body payload))
        data-dir (str (clj-str/replace (rseqc/correct-filepath filepath) #"/$" "") "/")
        log-path (fs-lib/join-paths workdir "log")
        response {:report (make-remote-link (format "%s/multiqc_report.html" workdir))
                  :log (make-remote-link log-path)}
        task-id (create-task! {:id             uuid
                               :name           name
                               :description    description
                               :payload        payload
                               :owner          owner
                               :plugin-name    v/plugin-name
                               :plugin-type    "ReportPlugin"
                               :plugin-version (:plugin-version plugin-context)
                               :response       response})
        result-dir (fs-lib/join-paths workdir "results")]
    (fs-lib/create-directories! result-dir)
    (spit log-path (json/write-str {:status "Running"
                                    :msg ""}))
    (update-process! task-id 0)
    (publish-event! "quartet_rseqc_report"
                    {:data-dir data-dir
                     :metadata metadata
                     :dest-dir workdir
                     :task-id task-id
                     :parameters {:name name
                                  :description description
                                  :plugin-name v/plugin-name
                                  :plutin-type "ReportPlugin"
                                  :plugin-version (:plugin-version plugin-context)}})
    response))

(defn- filter-mkdir-copy
  [fmc-datadir fmc-patterns fmc-destdir fmc-newdir]
  (let [files-keep (rseqc/batch-filter-files fmc-datadir fmc-patterns)
        files-keep-dir (fs-lib/join-paths fmc-destdir fmc-newdir)]
    (fs-lib/create-directories! files-keep-dir)
    (if (empty? files-keep)
      (log/warn (format "Cannot find any files with pattern %s, please check your data." fmc-patterns))
      (rseqc/copy-files! files-keep files-keep-dir {:replace-existing true}))))

(defn copy-files-to-dir
  [data-dir dest-dir]
  (filter-mkdir-copy (format "%s%s" data-dir "call-ballgown") [".*.txt"] dest-dir "ballgown")
  (filter-mkdir-copy (format "%s%s" data-dir "call-count") [".*gene_count_matrix.csv"] dest-dir "count")
  (filter-mkdir-copy (format "%s%s" data-dir "call-qualimapBAMqc") [".*tar.gz"] dest-dir "results/post_alignment_qc/bam_qc")
  (filter-mkdir-copy (format "%s%s" data-dir "call-qualimapRNAseq") [".*tar.gz"] dest-dir "results/post_alignment_qc/rnaseq_qc")
  (filter-mkdir-copy (format "%s%s" data-dir "call-fastqc") [".*.zip"] dest-dir "results/rawqc/fastqc")
  (filter-mkdir-copy (format "%s%s" data-dir "call-fastqscreen") [".*.txt"] dest-dir "results/rawqc/fastq_screen"))

(defn make-report!
  "Chaining Pipeline: filter-files -> copy-files -> merge_exp_file -> exp2qcdt -> multiqc."
  [{:keys [data-dir parameters metadata dest-dir task-id]}]
  (log/info "Generate quartet rnaseq report: " data-dir parameters metadata dest-dir)
  (let [metadata-file (fs-lib/join-paths dest-dir
                                         "results"
                                         "metadata.csv")
        parameters-file (fs-lib/join-paths dest-dir
                                           "results"
                                           "general-info.json")
        ballgown-dir (fs-lib/join-paths dest-dir "ballgown")
        count-dir (fs-lib/join-paths dest-dir "count")
        exp-fpkm-filepath (fs-lib/join-paths dest-dir "fpkm.csv")
        exp-count-filepath (fs-lib/join-paths dest-dir "count.csv")
        result-dir (fs-lib/join-paths dest-dir "results")
        log-path (fs-lib/join-paths dest-dir "log")
        config-path (fs-lib/join-paths dest-dir "results", "quartet_rnaseq_report.yaml")
        subdirs (rseqc/list-dirs data-dir)]
    (log/info "List subdirs: " subdirs)
    (try
      (fs-lib/create-directories! result-dir)
      (doseq [subdir subdirs]
        (copy-files-to-dir subdir dest-dir))
      (log/info "Merge gene experiment files from ballgown directory to a experiment table: " ballgown-dir exp-fpkm-filepath)
      (log/info "Merge gene experiment files from count directory to a experiment table: " count-dir exp-count-filepath)
      (update-process! task-id 10)
      (rseqc/merge-exp-files! (rseqc/list-files ballgown-dir {:mode "file"}) exp-fpkm-filepath)
      (rseqc/merge-exp-files! (rseqc/list-files count-dir {:mode "file"}) exp-count-filepath)
      (rseqc/write-csv! metadata-file metadata)
      ;;(decompression-tar files-qualimap-bam)
      ;;(decompression-tar files-qualimap-RNA)
      (doseq [files-qualimap-bam-tar (rseqc/batch-filter-files (fs-lib/join-paths dest-dir "results/post_alignment_qc/bam_qc") [".*tar.gz"])]
        (rseqc/decompression-tar files-qualimap-bam-tar))
      (doseq [files-qualimap-RNA-tar (rseqc/batch-filter-files (fs-lib/join-paths dest-dir "results/post_alignment_qc/rnaseq_qc") [".*tar.gz"])]
        (rseqc/decompression-tar files-qualimap-RNA-tar))
      (update-process! task-id 50)
      (rseqc/gen-multiqc-config config-path)
      (let [results (util/chain-fn-coll [(fn []
                                           (update-process! task-id 60)
                                           (rseqc/call-exp2qcdt! exp-fpkm-filepath exp-count-filepath metadata-file result-dir))
                                         (fn []
                                           (update-process! task-id 70)
                                           (spit parameters-file (json/write-str {"Report Name" (:name parameters)
                                                                                  "Description" (:description parameters)
                                                                                  "Report Tool" (format "%s-%s"
                                                                                                        (:plugin-name parameters)
                                                                                                        (:plugin-version parameters))
                                                                                  "Team" "Quartet Team"
                                                                                  "Date" (date)}))
                                           {:status "Success" :msg ""})
                                         (fn []
                                           (update-process! task-id 80)
                                           (rseqc/multiqc result-dir dest-dir {:config config-path
                                                                               :template "quartet_rnaseq_report"
                                                                               :title "Quartet RNA report"
                                                                               :env {:PATH (add-env-to-path "quartet-rseqc-report")}}))]
                                        (fn [result] (= (:status result) "Success")))
            status (:status (last results))
            msg (apply str (map :msg results))
            process (if (= status "Success") 100 -1)]
        (log/info (format "Running batch command: %s" (pr-str results)))
        (update-log-process! log-path {:status status
                                       :msg msg}
                             task-id process))
      (catch Exception e
        (update-process! task-id -1)
        (let [log (json/write-str {:status "Error" :msg (.toString e)})]
          (log/info "Status: " log)
          (spit log-path log))))))

(def events-init
  "Automatically called during startup; start event listener for quartet_rseqc_report events."
  (make-events-init "quartet_rseqc_report" make-report!))
