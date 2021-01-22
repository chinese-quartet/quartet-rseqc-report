#' Make a result directory for exp2qcdt
#'
#' @param workdir Working directory
make_directories <- function(workdir) {
  dirs <- sapply(c("performance_assessment", "rawqc", "post_alignment_qc", "quantification_qc", "simplified_report"), function(subdir) {
    return(file.path(workdir, subdir))
  })

  sapply(dirs, function(dest) {
    dir.create(dest, showWarnings = FALSE)
  })
}

#' Convert expression table to qc data table.
#'
#' @param exp_table_file Expression table file
#' @param phenotype_file Phenotype file
#' @param result_dir A directory for result files
#' @importFrom data.table fread
#' @importFrom data.table as.data.table
#' @importFrom data.table setnames
#' @importFrom data.table ":="
#' @importFrom dplyr %>%
#' @importFrom reshape2 melt
#' @importFrom reshape2 dcast
#' @importFrom utils combn
#' @export

exp2qcdt <- function(exp_table_file, count_table_file, phenotype_file, result_dir) {
  ref_data_dir <- paste(system.file(package = "exp2qcdt"), "/data", sep = "")
  # Global variable
  # TODO: This is not a good choice, maybe have another solution
  ref_data <<- read_ref_data(ref_data_dir)
  dt_fpkm <- fread(exp_table_file)
  dt_counts <- fread(count_table_file)
  dt_meta <- fread(phenotype_file)

  # Prepare directories
  make_directories(result_dir)

  if(colnames(dt_fpkm)[1] != 'gene_id'){
    colnames(dt_fpkm)[1] = 'gene_id'
  }
  
  if(!all(colnames(dt_counts) == colnames(dt_fpkm))){
    colnames(dt_fpkm) <- colnames(dt_counts)
  } 
  
  # expression data type must be numeric
  change_cols <- colnames(dt_fpkm[, !'gene_id'])
  dt_fpkm[, (change_cols):= lapply(.SD, as.numeric), .SDcols = change_cols]
  dt_counts[, (change_cols):= lapply(.SD, as.numeric), .SDcols = change_cols]
  
  dt_fpkm_log <- data.table(apply(dt_fpkm[, !'gene_id'], 2, function(x)(log2(x + 0.01))))
  dt_fpkm_log[, gene_id := dt_fpkm$gene_id]

  if (length(unique(dt_meta$sample)) < 2 & all(table(dt_meta[, 'sample']) < 2) ) {
    stop('At least two types of samples are required to calculate SNR')
  } else {
    one_group_out_list <- get_one_group(dt_fpkm_log, dt_counts, dt_meta, result_dir)
    abs_cor_median <- one_group_out_list[[1]]
    pt_abs_median_cor <- one_group_out_list[[2]]
    make_performance_plot(dt_fpkm, dt_fpkm_log, dt_counts, dt_meta, result_dir, abs_cor_median, pt_abs_median_cor)
  }
}

