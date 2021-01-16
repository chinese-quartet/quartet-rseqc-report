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
#' @importFrom cowplot ggdraw
#' @importFrom cowplot insert_xaxis_grob
#' @importFrom cowplot insert_yaxis_grob
#' @importFrom ggthemes theme_few
#' @importFrom ggplot2 theme_classic

make_performance_plot <- function(dt_fpkm_log, dt_counts, dt_meta, result_dir) {
  # two group which two replicates are need 
  
  sample_type_list <- dt_meta[['sample']] %>% unique()
  
  ### obtain relative expression level correlation table ------------------ 
  compare_combn <- data.table(combn(sample_type_list, 2))
  dt_rel_cor <- lapply(1:dim(compare_combn)[2], function(x){
    compare_name <- paste(compare_combn[[x]][2], '/', compare_combn[[x]][1], sep = '')
    dt_count_compare <- dt_counts[, c('gene_id',
                                      dt_meta[c(compare_combn[[x]][1], compare_combn[[x]][2]), on = .(sample)][['library']]), with = F]
    gene_list_com <- dt_count_compare[['gene_id']][apply(dt_count_compare[, !'gene_id'], 1, function(x){all(x >= 3)})]
    group_list1 <- dt_meta[sample == compare_combn[[x]][1]][['group']]
    group_list2 <- dt_meta[sample == compare_combn[[x]][2]][['group']]
    ratio_group <- outer(group_list2, group_list1, paste, sep = 'to')
    ratio_cor_group <- data.table(combn(ratio_group, 2))
    
    
    mat_rel_cor_per <- lapply(1:dim(ratio_cor_group)[2], function(y){
      ratio_name_1 <- ratio_cor_group[[y]][1]
      sample_id_a1 <- dt_meta[group == strsplit(ratio_name_1, 'to')[[1]][1]][['library']]
      sample_id_a2 <- dt_meta[group == strsplit(ratio_name_1, 'to')[[1]][2]][['library']]
      ratio_value_a <- dt_fpkm_log[gene_list_com, on = .(gene_id)][[sample_id_a1]] - dt_fpkm_log[gene_list_com, on = .(gene_id)][[sample_id_a2]]
      ratio_name_2 <- ratio_cor_group[[y]][2]
      sample_id_b1 <- dt_meta[group == strsplit(ratio_name_2, 'to')[[1]][1]][['library']]
      sample_id_b2 <- dt_meta[group == strsplit(ratio_name_2, 'to')[[1]][2]][['library']]
      ratio_value_b <- dt_fpkm_log[gene_list_com, on = .(gene_id)][[sample_id_b1]] - dt_fpkm_log[gene_list_com, on = .(gene_id)][[sample_id_b2]]
      cor_value <- cor(ratio_value_a, ratio_value_b)
      
      return(list(
        compare_group = compare_name,
        comapre_id = paste(ratio_name_1, ratio_name_2, sep = 'vs'),
        cor_value = cor_value,
        gene_num = length(gene_list_com)
      ))
    }) %>% rbindlist()
    
    return(mat_rel_cor_per)
  }) %>% rbindlist()
  
  ### relative replicates correlation of median group ----------------------
  # scatter plot and data between two replicates of one sample
  dt_rel_cor_median <- dt_rel_cor[with(dt_rel_cor, which.min(abs(dt_rel_cor$cor_value - median(dt_rel_cor$cor_value))))]
  output_rel_rep_res <- function(dt_rel_cor, dt_fpkm_log, dt_counts, dt_meta){
    compare_name_median <- dt_rel_cor_median[[comapre_id]]
    # obtain median relative correlation group
    sample_r1 <- strsplit(compare_name_median, 'vs')[[1]][1]
    sample_r2 <- strsplit(compare_name_median, 'vs')[[1]][2]
    sample_r1_1 <- strsplit(sample_r1, 'to')[[1]][1]
    sample_r1_2 <- strsplit(sample_r1, 'to')[[1]][2]
    sample_r2_1 <- strsplit(sample_r2, 'to')[[1]][1]
    sample_r2_2 <- strsplit(sample_r2, 'to')[[1]][2]
    
    # filter gene based on counts >= 3
    dt_count_compare <- dt_counts[, c('gene_id', dt_meta[c(sample_r1_1, sample_r1_2, sample_r2_1, sample_r2_2), on = .(group)][['library']]), with = F]
    gene_list_com <- dt_count_compare[['gene_id']][apply(dt_count_compare[, !'gene_id'], 1, function(x){all(x >= 3)})]
    
    # obtain relative expression table
    sample_id_a1 <- dt_meta[group == sample_r1_1][['library']]
    sample_id_a2 <- dt_meta[group == sample_r1_2][['library']]
    ratio_value_a <- dt_fpkm_log[gene_list_com, on = .(gene_id)][[sample_id_a1]] - dt_fpkm_log[gene_list_com, on = .(gene_id)][[sample_id_a2]]
    ratio_name_1 <- paste(sample_r1_1, '/', sample_r1_2, sep = '')
    sample_id_b1 <- dt_meta[group == sample_r2_1][['library']]
    sample_id_b2 <- dt_meta[group == sample_r2_2][['library']]
    ratio_value_b <- dt_fpkm_log[gene_list_com, on = .(gene_id)][[sample_id_b1]] - dt_fpkm_log[gene_list_com, on = .(gene_id)][[sample_id_b2]]
    ratio_name_2 <- paste(sample_r2_1, '/', sample_r2_2, sep = '')
    dt_median_rel_exp <- data.table(cbind(ratio_value_a, ratio_value_b))
    dt_median_rel_exp[, gene_id := gene_list_com]
    colnames(dt_median_rel_exp) <- c(ratio_name_1, ratio_name_2, 'gene_id')
    # output table 
    return(dt_median_rel_exp)
    }
  
  # relative expression table of two sample   
  dt_rel_scatter <- output_rel_rep_res(dt_rel_cor, dt_fpkm_log, dt_counts, dt_meta, result_dir)
  fwrite(dt_rel_scatter, file = paste(result_dir, "/performance_assessment/rel_median_sample_two_replicate_correlation.txt", sep = ""), sep = "\t" )
  
  ## output figure
  # relative correlation plot 
  xlab_rel_cor <- colnames(dt_rel_scatter)[1]
  ylab_rel_cor <- colnames(dt_rel_scatter)[2]
  colnames(dt_rel_scatter) <- c('replicate1', 'replicate2', 'gene_id')
  rel_cor_pt <- round(cor(dt_rel_scatter$replicate1, dt_rel_scatter$replicate2), digits = 3)
  pt_rel_median_cor <- ggplot2::ggplot(dt_rel_scatter, aes(x = replicate1, y = replicate2)) +
    geom_point(alpha = 0.8, size = 0.3, col = '#7BC8A4') +
    scale_fill_viridis_c(name = "density") +
    make_theme () +
    labs(title = 'Two Group: Relative Correlation',
         subtitle = paste('Correlation = ', cor_value_pt, ' N = (', gene_num, ')'),
         x = xlab_rel_cor,
         y = ylab_rel_cor)
  
  # mean relative correlation 
  logfc_corr_with_ref <- lapply(1:dim(compare_combn)[2], function(x){
    compare_name <- paste(compare_combn[[x]][2], '/', compare_combn[[x]][1], sep = '')
    corr_ref_compare <- corr_ref[compare_name, on = .(compare)]
    
    # get comapre column
    dt_count_compare <- dt_counts[, c('gene_id',
                                      dt_meta[c(compare_combn[[x]][1], compare_combn[[x]][2]), on = .(sample)][['library']]), with = F]
    gene_list_com <- dt_count_compare[['gene_id']][apply(dt_count_compare[, !'gene_id'], 1, function(x){all(x > 3)})]
    dt_count_compare_com <- dt_count_compare[gene_list_com, on = .(gene_id)]
    group_compare <- dt_meta[c(compare_combn[[x]][1], compare_combn[[x]][2]), on = .(sample)][['sample']]
    
    # DEG analysis was perform with edger pakcages
    deg_output <- DEGanalysis(dt_count_compare_com[, !'gene_id'], group_compare)
    deg_output[, gene_id := dt_count_compare_com[as.numeric(deg_output$gene), 'gene_id']]
    deg_output_ref <- deg_output[corr_ref_compare, on=.(gene_id==gene), nomatch = 0]
    
    # logfc correlation value
    logfc_cor_ref <- cor(deg_output_ref$logFC, deg_output_ref$meanlogFC)
    return(list(
      Batch = 'QC_test',
      DataQual = 'Test',
      compare = compare_name,
      corr_ref = logfc_cor_ref, 
      N_ref = dim(deg_output_ref)[1]
    ))
  }) %>% rbindlist()
  
  ### two group scatter plot figure and data
  # count mean relative correlation 
  dt_rel_cor_mean <- lapply(unique(dt_rel_cor$compare_group), function(x){
    mean_cor_value <- mean(dt_rel_cor[x, on = .(compare_group)][['cor_value']])
    return(list(
      compare = unique(dt_rel_cor[x, on = .(compare_group)][['compare_group']]),
      corr_FC = mean_cor_value,
      N_FC = unique(dt_rel_cor[x, on = .(compare_group)][['gene_num']])
    ))
  }) %>% rbindlist()
  
  # combine new relative correlation value and logfc correlation value from reference to reference data
  dt_cor_logfc_new <- logfc_corr_with_ref[dt_rel_cor_mean, on = .(compare)]
  dt_cor_logfc_combine <- rbind(refqc_202011_forplot, dt_cor_logfc_new)
  
  ### relative performance -----------------------------------
  # scaltter plot between logfc correlation based on reference logfc and relative expression correlation
  # relative performance correlation plot 
  pt_rel_fc <- ggplot2::ggplot(dt_cor_logfc_combine, aes(x = corr_ref, y = corr_FC)) +
    geom_point(alpha = 0.8, size = 0.3, col = 'steelblue') +
    #make_theme () +
    labs(title = 'Performance evaluation in relative expression',
         x = "Reference datasets based relative correlation",
         y = "Relative correlation")

  # relative performance output data
  fwrite(dt_cor_logfc_combine, file = paste(result_dir, "/performance_assessment/performance_evaluation_in_relative_expression.txt", sep = ""), sep = "\t")
  
  ### SNR performance -----------------------------------------
  ## obtain SNR results
  output_snr_res <- function(dt_fpkm, dt_counts, dt_meta){
    gene_list_snr <-  dt_counts[['gene_id']][apply(dt_counts[, !'gene_id'], 1, function(x){all(x >= 3)})]
    exp_design = (dt_meta[, .(library, group = sample)] %>% setkey(., library))
    dt_fpkm_f <- dt_fpkm[gene_list_snr, on = .(gene_id)]
    dt_fpkm_zscore <- data.table(t(apply(dt_fpkm_f[, !'gene_id'], 1, function(x){(x - mean(x))/sd(x)})))
    colnames(dt_fpkm_zscore) <- colnames(dt_fpkm_f[, !'gene_id'])
    pca_list <- get_pca_list(dt_fpkm_zscore, exp_design, dt_meta)
    return(pca_list)
  }
  
  dt_snr <- output_snr_res(dt_fpkm, dt_counts, dt_meta)
  snr_value <- dt_snr$SNR[1]
  snr_gene_num <- dt_snr$gene_num[1]
  
  ## figure of pca with snr
  pt_snr <- ggplot(dt_snr, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample), size = 2.5) +
    theme_few() +
    make_theme() +
    guides(shape = guide_legend(ncol = 1), color = guide_legend(ncol = 1, title.position = "top")) +
    theme(axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14)) +
    scale_fill_manual(values = c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745")) +
    scale_color_manual(values = c("#2f5c85", "#7BC8A4", "#FFC65D", "#F16745")) +
    ggtitle(paste("SNR: ", dt_snr$SNR[1], sep = "")) +
    labs(x = paste("PC1 (", dt_snr$PC1_ratio, "%)", sep = ""),
         y = paste("PC1 (", dt_snr$PC2_ratio, "%)", sep = ""))
  
  ## output snr table
  fwrite(dt_snr, file = paste(result_dir, "/performance_assessment/studydesign_performance_summary_more.txt", sep = ""), sep = "\t")
  
  ## abs expression correlation and snr
  rel_cor_median <- dt_rel_cor_median[['cor_value']]
  dt_snr_abs_rel_cor_new <- data.table(batch = 'QC_test', SNR = snr_value, LIR = abs_cor_median, 
                                       LRR2 = rel_cor_median, DataQual = 'Test')
  dt_snr_abs_rel_cor_combine <- rbind(qcintra_forplot, dt_snr_abs_rel_cor_new)
  fwrite(dt_snr_abs_rel_cor_combine, file = paste(result_dir, "/performance_assessment/snr_abs_rel_cor.txt", sep = ""), sep = "\t")
  
  ### output report figure-----------------------------
  pdf(file = paste(result_dir, "/simplified_report/", "figure2", ".pdf", sep = ""), 10, 3)
  plot_grid(pt_snr, pt_rel_median_cor, pt_abs_median_cor, 
            byrow = TRUE, ncol = 3)
  dev.off()
  
  #' Plot scatter and box combined plot
  plot_scatter_box <- function(dt_sb, var_x, var_y, col_g){
    pmain <- ggplot(dt_sb, aes(x = {{var_x}}, y = {{var_y}}, color = {{col_g}})) +
      geom_point()+
      scale_color_manual(values=c('#2f5c85','#7ba1c7','red')) +
      theme_few() +
      theme(legend.position = "none")
    
    xplot <- ggplot(dt_sb, aes(x = {{var_y}}, y = {{col_g}}, colour = {{col_g}})) + 
      geom_boxplot() +
      scale_color_manual(values=c('#2f5c85','#7ba1c7','red'))+
      coord_flip() +
      theme_classic()
    yplot <- ggplot(dt_sb, aes(x = {{col_g}}, y = {{var_x}}, colour = {{col_g}})) + 
      geom_boxplot() +
      scale_color_manual(values=c('#2f5c85','#7ba1c7','red'))+
      theme_classic()
    
    p1 <- insert_xaxis_grob(pmain, xplot, grid::unit(.2, "null"), position = "top")
    p2 <- insert_yaxis_grob(p1, yplot, grid::unit(.2, "null"), position = "right")
    pt_sb <- ggdraw(p2)
    return(pt_sb)
  }
  
  pt_snr_lir <- plot_scatter_box(dt_snr_abs_rel_cor_combine, var_x = 'SNR', var_y = 'LIR', col_g = 'DataQual')
  
  ### output summary table
  snr_rank <- c(rank(dt_snr_abs_rel_cor_combine$SNR)[which(dt_snr_abs_rel_cor_combine$batch == "QC_test")])
  abs_cor_rank <- c(rank(dt_snr_abs_rel_cor_combine$LIR)[which(dt_snr_abs_rel_cor_combine$batch == "QC_test")])
  rel_cor_rank <- c(rank(dt_snr_abs_rel_cor_combine$LRR2)[which(dt_snr_abs_rel_cor_combine$batch == "QC_test")])
  dt_metric_summary <- data.table(qc_metrics = c('Signal-to-Noise Ratio (SNR)', 'Relative correlation', 'Absolute correlation'),
                                  category = c('More groups', 'Two groups', 'One group'),
                                  value = c(snr_value, rel_cor_median, abs_cor_median),
                                  historical_value = c('14.45 ± 9.58', '0.493 ± 0.111', '0.973 ± 0.015'),
                                  rank = c(snr_rank, abs_cor_rank, rel_cor_rank))
  fwrite(dt_metric_summary, file = paste(result_dir, "/performance_assessment/qc_metric_summary.txt", sep = ""), sep = "\t")
  }