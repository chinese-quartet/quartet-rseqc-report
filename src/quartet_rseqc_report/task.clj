(ns quartet-rseqc-report.task
  (:require [quartet-rseqc-report.rseqc :as rseqc]
            [local-fs.core :as fs-lib]
            [tservice-core.plugins.env :refer [get-workdir add-env-to-path create-task! update-task!]]
            [tservice-core.plugins.util :as util]
            [clojure.data.json :as json]
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
  [{{:keys [name data_dir metadata_file description owner plugin-context]
     :or {description (format "Quality control report for %s" name)}
     :as payload} :body}]
  (log/info (format "Create a report %s with %s" name payload))
  (let [workdir (get-workdir)
        uuid (fs-lib/basename workdir)
        payload (merge {:description description} payload)
        data-dir (rseqc/correct-filepath data_dir)
        metadata-file (rseqc/correct-filepath metadata_file)
        log-path (fs-lib/join-paths workdir "log")
        response {:report (format "%s/multiqc_report.html" workdir)
                  :log log-path}
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
                     :metadata-file metadata-file
                     :dest-dir workdir
                     :task-id task-id
                     :metadata {:name name
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

(defn make-report!
  "Chaining Pipeline: filter-files -> copy-files -> merge_exp_file -> exp2qcdt -> multiqc."
  [{:keys [datadir parameters metadata-file dest-dir task-id]}]
  (log/info "Generate quartet rnaseq report: " datadir parameters metadata-file dest-dir)
  (let [parameters-file (fs-lib/join-paths dest-dir
                                           "results"
                                           "general-info.json")
        ballgown-dir (fs-lib/join-paths dest-dir "ballgown")
        count-dir (fs-lib/join-paths dest-dir "count")
        exp-fpkm-filepath (fs-lib/join-paths dest-dir "fpkm.csv")
        exp-count-filepath (fs-lib/join-paths dest-dir "count.csv")
        result-dir (fs-lib/join-paths dest-dir "results")
        log-path (fs-lib/join-paths dest-dir "log")
        config-path (fs-lib/join-paths dest-dir "results", "quartet_rnaseq_report.yaml")]
    (try
      (fs-lib/create-directories! result-dir)
      (log/info "Merge gene experiment files from ballgown directory to a experiment table: " ballgown-dir exp-fpkm-filepath)
      (log/info "Merge gene experiment files from count directory to a experiment table: " count-dir exp-count-filepath)
      (filter-mkdir-copy datadir [".*ballgown/.*.txt"] dest-dir "ballgown")
      (filter-mkdir-copy datadir [".*count/.*gene_count_matrix.csv"] dest-dir "count")
      (filter-mkdir-copy datadir [".*qualimapBAMqc/.*tar.gz"] dest-dir "results/post_alignment_qc/bam_qc")
      (filter-mkdir-copy datadir [".*qualimapRNAseq/.*tar.gz"] dest-dir "results/post_alignment_qc/rnaseq_qc")
      (filter-mkdir-copy datadir [".*fastqc/.*.zip"] dest-dir "results/rawqc/fastqc")
      (filter-mkdir-copy datadir [".*fastqscreen/.*.txt"] dest-dir "results/rawqc/fastq_screen")
      (update-process! task-id 10)
      (rseqc/merge-exp-files! (rseqc/list-files ballgown-dir {:mode "file"}) exp-fpkm-filepath)
      (rseqc/merge-exp-files! (rseqc/list-files count-dir {:mode "file"}) exp-count-filepath)
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
                                                                               :env {:PATH (add-env-to-path "quartet-rnaseq-report")}}))]
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
