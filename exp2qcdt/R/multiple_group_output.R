#' Make performance figure
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 scale_fill_viridis_c
#' @importFrom cowplot ggdraw
#' @importFrom cowplot insert_xaxis_grob
#' @importFrom cowplot insert_yaxis_grob
#' @importFrom cowplot plot_grid
#' @importFrom ggthemes theme_few
#' @importFrom ggplot2 theme_classic
#' @importFrom dplyr bind_rows
#' @importFrom scales rescale
#' @importFrom data.table rbindlist
#' @importFrom data.table fwrite
#' @importFrom data.table ":="
#' @importFrom data.table "setDF"
#' @importFrom data.table "setDT"
#' @importFrom utils combn

make_performance_plot <- function(dt_fpkm, dt_fpkm_log, dt_counts, dt_meta, result_dir, 
                                  abs_cor_median, pt_abs_median_cor) {
  
  # import reference data
  dt_ref_qc_metrics_value <- ref_data$ref_qc_metrics_value
  dt_ref_fc_value <- ref_data$ref_fc_value
  
  # two group which two replicates are need 
  sample_type_list <- dt_meta[['sample']] %>% unique()
  
  ### D5/D6, F7/D6, M8/D6 log2FC correlation with reference data (dt_ref_fc_value)  ------------------ 
  compare_combn <- data.table(combn(sample_type_list, 2))
  
  # test data logfc 
  dt_fc_test <- do.call(rbind, lapply(list(c('D5', 'D6'), c('F7', 'D6'), c('M8', 'D6')), function(x){
    compare_name <- paste(x[1], '/', x[2], sep = '')
    
    ### at least tow replicate counts >= 3 
    dt_detect_gene <-  data.table(apply(dt_counts[, dt_meta[x[1], on = .(sample)][['library']], with = F], 1, function(x){length(which(x >= 3)) >= 2}),
                                  apply(dt_counts[, dt_meta[x[2], on = .(sample)][['library']], with = F], 1, function(x){length(which(x >= 3)) >= 2}))
    gene_list_com <- dt_counts[['gene_id']][apply(dt_detect_gene, 1, function(x){all(x)})]
    
    dt_count_compare <- dt_counts[, c('gene_id',
                                      dt_meta[c(x[1], x[2]), on = .(sample)][['library']]), with = F]
    dt_count_compare_com <- dt_count_compare[gene_list_com, on = .(gene_id)]
    group_compare <- dt_meta[c(x[1], x[2]), on = .(sample)][['sample']]
    group_compare[group_compare == x[1]] <- 'ZZZ'
    group_compare[group_compare == x[2]] <- 'AAA'
    
    # DEG analysis was perform with edger pakcages
    deg_output <- DEGanalysis(dt_count_compare_com[, !'gene_id'], group_compare)
    deg_output[, gene_id := dt_count_compare_com[as.numeric(deg_output$gene), 'gene_id']]
    deg_output[, compare := compare_name]
    deg_output_d <- deg_output[,c('gene_id', 'compare', 'logFC'), with = FALSE]
    colnames(deg_output_d) <- c('gene', 'compare', 'meanlogFC')
    return(deg_output_d)
  }))
  
  # get reference data set compared group
  dt_fc_test[, gene_compare := paste(gene, compare, sep = '_')]
  dt_ref_fc_value[, gene_compare := paste(gene, compare, sep = '_')]
  dt_ref_fc_test <- dt_fc_test[dt_ref_fc_value, on = 'gene_compare', nomatch = 0]
  dt_ref_fc_test_d <- dt_ref_fc_test[, c('gene', 'compare', 'meanlogFC', 'i.meanlogFC'), with = FALSE]
  colnames(dt_ref_fc_test_d) <- c('gene', 'compare', 'meanlogFC_test', 'meanlogFC_ref')
  
  # log2fc correlation output data
  cor_log2fc <- format(round(cor(dt_ref_fc_test_d$meanlogFC_test, dt_ref_fc_test_d$meanlogFC_ref), digits = 3), nsmall = 3)
  dt_ref_fc_test_d[, meanlogFC_test := round(meanlogFC_test, digits = 3)]
  dt_ref_fc_test_d[, cor := cor_log2fc][, gene_num := dim(dt_ref_fc_test_d)[1]]
  fwrite(dt_ref_fc_test_d, file = paste(result_dir, "/performance_assessment/logfc_cor_ref_test.txt", sep = ""), sep = "\t")
  
  # log2fc correlation output figure
  pt_logfc_cor <- ggplot2::ggplot(dt_ref_fc_test_d, aes(x = meanlogFC_ref, y = meanlogFC_test, color = compare)) +
    geom_point(alpha = 0.8, size = 0.3) +
    theme_few() + 
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    scale_fill_viridis_c(name = "density") +
    labs(title = 'LogFC Correlation',
         subtitle = paste('Correlation: ', cor_log2fc, ' (N = ', dim(dt_ref_fc_test_d)[1], ')', sep = ''),
         x = 'meanlogFC (Reference)',
         y = 'meanlogFC (Test)')
  
  ### SNR performance -----------------------------------------
  ## obtain SNR results
  output_snr_res <- function(dt_fpkm_log, dt_counts, dt_meta){
    dt_detect_gene <- do.call(cbind, lapply(unique(dt_meta$sample), function(x){
      detect_res <- apply(dt_counts[, dt_meta[sample == x][['library']], with = F], 1, function(x){length(which(x >= 3)) >= 2})
      return(
        detect_res = detect_res
        ) 
    }))
    
    gene_list_snr <- dt_counts[['gene_id']][apply(dt_detect_gene, 1, function(x){any(x)})]
    exp_design = (dt_meta[, .(library, group = sample)] %>% setkey(., library))
    dt_fpkm_f <- dt_fpkm_log[gene_list_snr, on = .(gene_id)]
    dt_fpkm_zscore <- data.table(t(apply(dt_fpkm_f[, dt_meta$library, with = F], 1, function(x){(x - mean(x))/sd(x)})))
    pca_list <- get_pca_list(dt_fpkm_zscore, exp_design, dt_meta)
    return(pca_list)
  }
  
  dt_snr <- output_snr_res(dt_fpkm_log, dt_counts, dt_meta)
  snr_gene_num <- dt_snr$gene_num[1]
  
  ## figure of pca with snr
  pt_snr <- ggplot(dt_snr, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample), size = 2.5, show.legend = FALSE) +
    theme_few() +
    guides(shape = guide_legend(ncol = 1), color = guide_legend(ncol = 1, title.position = "top")) +
    scale_fill_manual(values = c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745")) +
    scale_color_manual(values = c("#2f5c85", "#7BC8A4", "#FFC65D", "#F16745")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(
      title = paste("SNR: ", dt_snr$SNR[1], ' (N = ', dt_snr$gene_num[1], ')', sep = ""),
      x = paste("PC1 (", dt_snr$PC1_ratio, "%)", sep = ""),
      y = paste("PC1 (", dt_snr$PC2_ratio, "%)", sep = ""))
  
  ## output snr table
  dt_snr$PC1 <- round(dt_snr$PC1 , digits = 3)
  dt_snr$PC2 <- round(dt_snr$PC2 , digits = 3)
  fwrite(dt_snr, file = paste(result_dir, "/performance_assessment/pca_with_snr.txt", sep = ""), sep = "\t")
  
  ### output report data and figure-----------------------------
  data.table::setDF(dt_ref_qc_metrics_value)
  test_metrics_value <- c('QC_test', as.character(dt_snr$SNR[1]), cor_log2fc, rep(NA, 6))
  dt_ref_qc_metrics_value[nrow(dt_ref_qc_metrics_value) + 1, ] <- test_metrics_value 
  data.table::setDT(dt_ref_qc_metrics_value)
  dt_ref_qc_metrics_value[, SNR := as.numeric(SNR)][, RC := as.numeric(RC)][batch == 'QC_test', total_score := format(round(sqrt(SNR*RC), digits = 3), nsmall = 3)]
  dt_ref_qc_metrics_value[batch == 'QC_test', group := 'Query'][batch != 'QC_test', group := 'Reference']
  dt_ref_qc_metrics_value[, total_score := as.numeric(total_score)]
  
  # log2fc correlation value and snr scatter plot
  pt_snr_rc_cor <- plot_scatter_box(dt_ref_qc_metrics_value, var_x = 'SNR', var_y = 'RC', 
                                     col_g = 'group', xlab = 'SNR', ylab = 'logfc correlation', 
                                     title_lab = 'Performance evaluation')

  pdf(file = paste(result_dir, "/simplified_report/", "figure1", ".pdf", sep = ""), 12, 4)
  pt_fig1 <- plot_grid(pt_logfc_cor, pt_snr, pt_snr_rc_cor, byrow = TRUE, ncol = 3)
  print(pt_fig1)
  dev.off()
  
  ### output summary table ---
  rank_len = dim(dt_ref_qc_metrics_value)[1]
  snr_val <- dt_ref_qc_metrics_value$SNR
  names(snr_val) <- dt_ref_qc_metrics_value$batch
  snr_rank <-  which(names(sort(snr_val, decreasing = TRUE)) == 'QC_test')
  rc_cor_val <- dt_ref_qc_metrics_value$RC
  names(rc_cor_val) <- dt_ref_qc_metrics_value$batch
  rc_cor_rank <-  which(names(sort(rc_cor_val, decreasing = TRUE)) == 'QC_test')
  total_score_val <- dt_ref_qc_metrics_value$total_score
  names(total_score_val) <- dt_ref_qc_metrics_value$batch
  total_score_rank <-  which(names(sort(total_score_val, decreasing = TRUE)) == 'QC_test')
  
  dt_metric_summary <- data.table(
    qc_metrics = c('Signal-to-Noise Ratio (SNR)', 'Relative Correlation with Reference Datasets (RC) ', 'Total Score'),
    value = c(as.character(dt_snr$SNR[1]), cor_log2fc, as.numeric(dt_ref_qc_metrics_value[batch == 'QC_test'][['total_score']])),
    historical_value = c('19.505 ± 7.039', '0.950 ± 0.028', '4.238 ± 0.849'),
    rank = c(paste(snr_rank, '/', rank_len, sep = ''), paste(rc_cor_rank, '/', rank_len, sep = ''), paste(total_score_rank, '/', rank_len, sep = '')))
  
  colnames(dt_metric_summary) <- c('qc_metrics', 'value', 'historical_value', 'rank')
  
  # qc metrics summary and rank
  fwrite(dt_metric_summary, file = paste(result_dir, "/performance_assessment/qc_metrics_summary.txt", sep = ""), sep = "\t")
  
  # rank qc metrics value and output
  dt_ref_qc_metrics_value_s <- dt_ref_qc_metrics_value[order(dt_ref_qc_metrics_value$total_score, decreasing = TRUE)]
  dt_ref_qc_metrics_value_s[, rank := 1: dim(dt_ref_qc_metrics_value_s)[1]]
  dt_ref_qc_metrics_value_s[rank < dim(dt_ref_qc_metrics_value_s)[1]/5, performance := 'Great']
  dt_ref_qc_metrics_value_s[dim(dt_ref_qc_metrics_value_s)[1]/5 <= rank & rank <= dim(dt_ref_qc_metrics_value_s)[1]*1/2, performance := 'Good']
  dt_ref_qc_metrics_value_s[dim(dt_ref_qc_metrics_value_s)[1]/5 < rev(rank) & rev(rank) <= dim(dt_ref_qc_metrics_value_s)[1]*1/2, performance := 'Fair']
  dt_ref_qc_metrics_value_s[rev(rank) < dim(dt_ref_qc_metrics_value_s)[1]/5, performance := 'Bad']
  fwrite(dt_ref_qc_metrics_value_s, paste(result_dir, "/performance_assessment/quality_score.txt", sep = ""), sep = "\t")
  
  ### quality score plot ---
  make_score_figure(result_dir, dt_ref_qc_metrics_value)
  }
