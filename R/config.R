#' Read reference data
#'
#' @param ref_data_dir A root directory for reference data
#' @return The reference data as a list.
#' @importFrom data.table fread
#' @export
read_ref_data <- function(ref_data_dir) {
  # Read reference data
  #
  # Args:
  #   None
  #
  # Returns:
  #   ref_data
  ref_data <- list()
  ref_data$highqual_batch <- c("R_ILM_L5_B2", "R_BGI_L7_B1", "R_ILM_L6_B1", "P_ILM_L6_B1", "R_ILM_L1_B1",
                      "R_BGI_L6_B1", "R_BGI_L3_B1", "R_ILM_L2_B2", "R_ILM_L4_B2", "R_ILM_L8_B1",
                      "R_ILM_L4_B3", "P_BGI_L6_B1", "P_BGI_L3_B1", "P_ILM_L1_B1", "P_ILM_L8_B1", "Test")
  ref_data$ref_meta <- fread(paste0(ref_data_dir, "/ref_metadata.csv"))
  ref_data$ref_data <- fread(paste0(ref_data_dir, "/ref_detect_sets.csv"))
  ref_data$ref_detected_gene_per <- fread(paste0(ref_data_dir, "/performance_evaluation_gene_detection.csv"))
  ref_data$ref_range <- fread(paste0(ref_data_dir, "/ref_reletive_exp_range.csv"))
  ref_data$snr_lst_gene <- fread(paste0(ref_data_dir, "/ref_snr_lst_gene.csv"))[["lst.gene"]]
  ref_data$ref_rel_exp_per <- fread(paste0(ref_data_dir, "/performance_evaluation_reletive_expression.csv"))
  ref_data$SD_ref <- fread(paste0(ref_data_dir, "/ref_studydesign_qc_summary.csv"))
  ref_data$ref_degs <- fread(paste0(ref_data_dir, "/ref_degs_all.csv"))
  ref_data$ref_degs_per <- fread(paste0(ref_data_dir, "/performance_evaluation_DEG.csv"))

  return(ref_data)
}

#' This function is for identify DEG from a given matrix
#'
#' @param expr_mat should a matrix in stead of a data frame, with column for samples and rows for signature(gene/protein/metabolics...)
#'                 For transcriptome, a matrix after log transformation is recommended.
#' @param group    should be a factor whose length is identical to the number of the columns in expr_mat,
#'                 describing the group information of each column in expr_mat
#' @param thr_filter genes with average no larger than thr_filter in both groups will not be considered as DEGs
#' @return A statistics table for deg
#' @importFrom stats sd
#' @importFrom stats t.test
#' @importFrom stats p.adjust
#' @importFrom dplyr %>%
#' @importFrom data.table rbindlist
#' @export
DEG_analysis <- function(expr_mat, group, thr_filter = log2(0.5), thr_FC = 2, thr_p = 0.05) {
  group_A <- levels(group)[1]
  group_B <- levels(group)[2]
  col_A <- which(as.numeric(group) == 1)
  col_B <- which(as.numeric(group) == 2)
  expr_mat_filt <- expr_mat[apply(expr_mat, 1, sd) > 0,]


  dt_DEG <- apply(expr_mat_filt, 1, function(x) {

    mean <- mean(x[c(col_A, col_B)])
    sd <- sd(x[c(col_A, col_B)])
    mean_A <- mean(x[col_A])
    mean_B <- mean(x[col_B])
    logFC <- mean_B - mean_A

    if ((mean_A > thr_filter) | (mean_B > thr_filter)) {
      p.value <- t.test(x[col_A], x[col_B], var.equal = TRUE)$p.value
      if (p.value < thr_p) {
        DEGtype <- "not DEG"
        if (logFC > log2(thr_FC)) DEGtype <- "Up-regulated"
        if (logFC < (-log2(thr_FC))) DEGtype <- "Down-regulated"
      } else {
        DEGtype <- "not DEG"
      }
    } else {
      p.value <- NA
      DEGtype <- "Low-expressed"
    }


    return(list(
      mean = mean,
      sd = sd,
      mean_A = mean_A,
      mean_B = mean_B,
      logFC = logFC,
      p.value = p.value,
      DEGtype = DEGtype
    ))

  }) %>% rbindlist()

  dt_DEG$group_A <- group_A
  dt_DEG$group_B <- group_B
  dt_DEG$gene <- rownames(expr_mat_filt)


  dt_DEG_passlowflit <- dt_DEG[DEGtype != "Low-expressed"]
  dt_DEG_passlowflit$p.value.adj <- p.adjust(dt_DEG_passlowflit$p.value, method = "fdr")
  dt_DEG.passLowfilt_sorted <- dt_DEG_passlowflit[order(ifelse(DEGtype %in% c("Up-regulated", "Down-regulated"), 0, 1),
                                                        - abs(logFC))]

  dt_DEG.lowExpr <- dt_DEG[DEGtype == "Low-expressed"]

  dt_DEG <- rbindlist(list(dt_DEG.passLowfilt_sorted, dt_DEG.lowExpr), use.names = TRUE, fill = TRUE)


  dt_DEG_foroutput <- dt_DEG[, .(group_A, group_B, gene, mean, mean_A, mean_B, sd, logFC, p.value, p.value.adj, DEGtype)]
  return(dt_DEG_foroutput)
}

#' Make a custom theme
#'
#' @return A custom theme
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_blank
#' @export
make_theme <- function() {
  custom_theme <- theme(plot.background = element_rect(colour = "white"),
                        axis.title.y = element_text(size = 16),
                        axis.title.x = element_text(size = 16),
                        axis.text.y = element_text(colour = "black", size = 16),
                        axis.text.x = element_text(colour = "black", size = 16),
                        title = element_text(colour = "black", size = 16),
                        panel.background = element_rect(fill = "white"),
                        panel.grid = element_blank(),
                        strip.text = element_text(size = 16))
  return(custom_theme)
}

#' Combine sd summary based on number of sample types
#'
#' @param result_dir A directory for result files
#' @param sample_num The number of samples
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @export
combine_sd_summary_table <- function(result_dir, sample_num) {
  if (sample_num > 2) {
    summary_one <- fread(paste(result_dir, "/performance_assessment/studydesign_performance_summary_one.txt", sep = ""))
    summary_more <- fread(paste(result_dir, "/performance_assessment/studydesign_performance_summary_more.txt", sep = ""))
    summary_table <- rbind(summary_one, summary_more)
    fwrite(summary_table, file = paste(result_dir, "/simplified_report/studydesign_performance_summary.txt", sep = ""), sep = "\t")
    fwrite(summary_table, file = paste(result_dir, "/simplified_report/SD_performace_table.txt", sep = ""), sep = "\t")
  } else {
    summary_one <- fread(paste(result_dir, "/performance_assessment/studydesign_performance_summary_one.txt", sep = ""))
    summary_table <- summary_one
    fwrite(summary_table, file = paste(result_dir, "/simplified_report/studydesign_performance_summary.txt", sep = ""), sep = "\t")
    fwrite(summary_table, file = paste(result_dir, "/simplified_report/SD_performace_table.txt", sep = ""), sep = "\t")
  }
}

#' Make the performance score figure
#'
#' @param result_dir A directory for result files
#' @param sample_num The number of samples
#' @importFrom dplyr %>%
#' @importFrom data.table data.table
#' @importFrom data.table fwrite
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 theme_void
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 arrow
#' @importFrom ggplot2 margin
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom grid unit
#' @export
make_score_figure <- function(result_dir, sample_num) {
  if (sample_num == 1) {
    score_table <- detected_gene_performance_mean[SD_performance_mean_one, on = "Batch"]
  } else if (sample_num == 2) {
    score_table <- detected_gene_performance_mean[rel_exp_performance_mean, on = "Batch"][degs_performance_mean, on = "Batch"][SD_performance_mean_one, on = "Batch"]
  } else if (sample_num > 2) {
    SD_performance_mean <- SD_performance_mean_one[SD_performance_mean_more, on = "Batch"]
    score_table <- detected_gene_performance_mean[rel_exp_performance_mean, on = "Batch"][degs_performance_mean, on = "Batch"][SD_performance_mean, on = "Batch"]
  }

  score_table_list <- colnames(score_table)[-1] %>% lapply(., function(x) {
    scale_value <- (score_table[[x]] - min(score_table[[x]])) / (max(score_table[[x]]) - min(score_table[[x]]))
    return({
      scale_value
    })
  })

  score_table_scale_mat <- data.table(do.call(cbind, score_table_list))
  score_table_scale_mat[, "Benchmarkingscore"] <- apply(score_table_scale_mat, 1, function(x) { mean(x) })
  score_table_scale_mat_m <- data.table(score_table$Batch, score_table_scale_mat)
  colnames(score_table_scale_mat_m) <- c(colnames(score_table), "Benchmarkingscore")
  score_table_scale_mat_m_order <- score_table_scale_mat_m[order(score_table_scale_mat_m$Benchmarkingscore, decreasing = TRUE),]
  fwrite(score_table_scale_mat_m_order, file = paste(result_dir, "/performance_assessment/performance_score.txt", sep = ""), sep = "\t")

  # S5-1 performance score figure
  dt_pscore <- data.table(cbind("score", score_table_scale_mat_m_order[, .(Benchmarkingscore, Batch)]))
  setnames(dt_pscore, "V1", "type")
  dt_pscore$Benchmarkingscore <- as.character(round(dt_pscore$Benchmarkingscore, digits = 2))
  test_score <- dt_pscore[.("Test"), on = .(Batch)][["Benchmarkingscore"]]
  pdf(paste(result_dir, "/simplified_report/performance_score.pdf", sep = ""), 4, 4)
  pt <- ggplot(dt_pscore, aes(x = Benchmarkingscore, y = type, fill = Benchmarkingscore)) +
    geom_tile(color = "white", show.legend = FALSE) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(9, "RdYlGn"))(16)) +
    annotate(geom = "curve", x = test_score,
             y = 2.5, xend = test_score, curvature = 0,
             yend = 1.5, arrow = arrow(angle = 45, length = unit(9, "mm"), type = "closed"), color = "grey") +
    annotate(geom = "text", x = dt_pscore[.("Test"), on = .(Batch)][["Benchmarkingscore"]],
             y = 2.2, label = test_score, hjust = "center", size = 10, fontface = "bold") +
    theme_void() +
    theme(plot.margin = margin(3, 0, 4, 0, "cm"))
  print(pt)
  dev.off()
}