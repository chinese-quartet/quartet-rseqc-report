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
#' @export
exp2qcdt <- function(exp_table_file, phenotype_file, result_dir) {
  ref_data_dir <- paste(system.file(package = "exp2qcdt"), "/data", sep = "")
  ref_data <- read_ref_data(ref_data_dir)

  dt_exp <- fread(exp_table_file)
  dt_meta <- fread(phenotype_file)

  # Prepare directories
  make_directories(result_dir)
  dt_exp_melt <- as.data.table(melt(dt_exp, id = 'GENE_ID'))
  setnames(dt_exp_melt, 1:3, c('gene', 'library', 'fpkm'))
  dt_exp_melt$gene <- as.character(dt_exp_melt$gene)
  dt_exp_melt$library <- as.character(dt_exp_melt$library)

  # get pairs table
  lst_library_all <- dt_meta$library
  dt_pairs <- combn(lst_library_all, 2) %>% t %>% as.data.table
  setnames(dt_pairs, 1:2, c('library_A', 'library_B'))
  dt_pairs$sampleA <- dt_meta[dt_pairs$library_A, on = 'library']$sample
  dt_pairs$sampleB <- dt_meta[dt_pairs$library_B, on = 'library']$sample
  dt_pairs$sample_type <- apply(dt_pairs, 1, function(x) {
    if (x['sampleA'] == x['sampleB']) return('Intra-sample')
    return('xross-sample')
  })
  dt_pairs$sample_type <- factor(dt_pairs$sample_type, levels = c('Intra-sample', 'xross-sample'))
  dt_exp_annot <- merge(dt_meta, dt_exp_melt, by = 'library')
  exp_fpkm <- dcast(dt_exp_annot, gene ~ library, value.var = 'fpkm') %>% data.frame(row.names = 1) %>% as.matrix()
  exp_fpkm_log <- log(exp_fpkm + 0.01)

  if (length(unique(dt_meta$sample)) < 1) {
    stop('There is no quantitative qc result')
  } else if (length(unique(dt_meta$sample)) == 1) {
    get_one_group(dt_exp_annot, dt_pairs, result_dir)
    combine_sd_summary_table(result_dir, sample_num = 1)
    make_score_figure(result_dir, sample_num = 1)
  } else if (length(unique(dt_meta$sample)) == 2) {
    get_one_group(dt_exp_annot, dt_pairs, result_dir)
    get_two_group(dt_exp_annot, dt_pairs, result_dir)
    combine_sd_summary_table(result_dir, sample_num = 2)
    make_score_figure(result_dir, sample_num = 2)
  } else if (length(unique(dt_meta$sample)) > 2) {
    get_one_group(dt_exp_annot, dt_pairs, result_dir)
    get_two_group(dt_exp_annot, dt_pairs, result_dir)
    get_more_group(exp_fpkm_log, dt_meta, result_dir)
    combine_sd_summary_table(result_dir, sample_num = 4)
    make_score_figure(result_dir, sample_num = 4)
  }
}

