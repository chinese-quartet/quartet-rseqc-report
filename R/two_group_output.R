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
                                      dt_meta[c(compare_combn[[x]][1], compare_combn[[x]][2]), on = .(sample)][['sample_id']]), with = F]
    gene_list_com <- dt_count_compare[['gene_id']][apply(dt_count_compare[, !'gene_id'], 1, function(x){all(x > 3)})]
    group_list1 <- dt_meta[sample == compare_combn[[x]][1]][['group']]
    group_list2 <- dt_meta[sample == compare_combn[[x]][2]][['group']]
    ratio_group <- outer(group_list2, group_list1, paste, sep = 'to')
    ratio_cor_group <- data.table(combn(ratio_group, 2))
    
    
    mat_rel_cor_per <- lapply(1:dim(ratio_cor_group)[2], function(y){
      ratio_name_1 <- ratio_cor_group[[y]][1]
      sample_id_a1 <- dt_meta[group == strsplit(ratio_name_1, 'to')[[1]][1]][['sample_id']]
      sample_id_a2 <- dt_meta[group == strsplit(ratio_name_1, 'to')[[1]][2]][['sample_id']]
      ratio_value_a <- dt_fpkm_log[gene_list_com, on = .(gene_id)][[sample_id_a1]] - dt_fpkm_log[gene_list_com, on = .(gene_id)][[sample_id_a2]]
      ratio_name_2 <- ratio_cor_group[[y]][2]
      sample_id_b1 <- dt_meta[group == strsplit(ratio_name_2, 'to')[[1]][1]][['sample_id']]
      sample_id_b2 <- dt_meta[group == strsplit(ratio_name_2, 'to')[[1]][2]][['sample_id']]
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
  
  ### relative replicates correlation ----------------------
  # scatter plot and data between two replicates of one sample
  compare_name_median <- dt_rel_cor[with(dt_rel_cor, which.min(abs(dt_rel_cor$cor_value - median(dt_rel_cor$cor_value))))][['comapre_id']]
  sample_r1 <- strsplit(compare_name_median, 'vs')[[1]][1]
  sample_r2 <- strsplit(compare_name_median, 'vs')[[1]][2]
  sample_r1_1 <- strsplit(sample_r1, 'to')[[1]][1]
  sample_r1_2 <- strsplit(sample_r1, 'to')[[1]][2]
  sample_r2_1 <- strsplit(sample_r2, 'to')[[1]][1]
  sample_r2_2 <- strsplit(sample_r2, 'to')[[1]][2]
  
  # filter gene based on counts >= 3
  dt_count_compare <- dt_counts[, c('gene_id', dt_meta[c(sample_r1_1, sample_r1_2, sample_r2_1, sample_r2_2), on = .(group)][['sample_id']]), with = F]
  gene_list_com <- dt_count_compare[['gene_id']][apply(dt_count_compare[, !'gene_id'], 1, function(x){all(x > 3)})]
  
  sample_id1 <- dt_meta[group == sample_1][['sample_id']]
  sample_id2 <- dt_meta[group == sample_2][['sample_id']]
  gene_list_1 <- dt_counts[dt_counts[[sample_id1]] >= 3][['gene_id']]
  gene_list_2 <- dt_counts[dt_counts[[sample_id2]] >= 3][['gene_id']]
  gene_list_com <- intersect(gene_list_1, gene_list_2)
  cor_vlaue <- round(dt_abs_cor[with(dt_abs_cor, which.min(abs(dt_abs_cor$abs_cor - median(dt_abs_cor$abs_cor))))][['abs_cor']], digits = 3)
  gene_num <- dt_abs_cor[with(dt_abs_cor, which.min(abs(dt_abs_cor$abs_cor - median(dt_abs_cor$abs_cor))))][['gene_num']]
  dt_one_group_out <- data.table(batch = 'QC_test', LIR = cor_vlaue)
  
  # output correlation plot and data
  dt_one_cor <- dt_fpkm_log[gene_list_com, on = .(gene_id)][, c('gene_id', sample_id1, sample_id2), with = F]
  fwrite(dt_one_cor, file = paste(result_dir, "/performance_assessment/median_sample_two_replicate_correlation.txt", sep = ""), sep = "\t")
  dt_scatter <- dt_one_cor
  colnames(dt_scatter) <- c('gene_id', 'replicate1', 'replicate2')
  
  # abs correlation plot 
  pt_abs_cor <- ggplot2::ggplot(dt_scatter, aes(x = replicate1, y = replicate2)) +
    geom_point(alpha = 0.8, size = 0.3, col = 'steelblue') +
    scale_fill_viridis_c(name = "density") +
    make_theme () +
    labs(title = 'One Group: Absolute Correlation',
         subtitle = paste('Correlation = ', cor_vlaue, ' N = (', gene_num, ')'),
         x = sample_1,
         y = sample_2)
  
  # output abs correlation plot 
  pdf(file = paste(result_dir, "/simplified_report/", "median_sample_two_replicate_correlation", ".pdf", sep = ""), 6, 5)
  print(pt)
  dev.off()
  
  ### obtain logfc correlation value between new data and reference data
  corr_ref <- fread("./data/TableS2_ReferenceDatasets.csv", drop = 'V1')
  refqc_202011_forplot <- readRDS("./data/refqc_202011_forplot.rds")
  
  logfc_corr_with_ref <- lapply(1:dim(compare_combn)[2], function(x){
    compare_name <- paste(compare_combn[[x]][2], '/', compare_combn[[x]][1], sep = '')
    corr_ref_compare <- corr_ref[compare_name, on = .(compare)]
    
    # get comapre column
    dt_count_compare <- dt_counts[, c('gene_id',
                                      dt_meta[c(compare_combn[[x]][1], compare_combn[[x]][2]), on = .(sample)][['sample_id']]), with = F]
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
  
  # output abs correlation plot 
  pdf(file = paste(result_dir, "/simplified_report/", "performance_evaluation_in_relative_expression", ".pdf", sep = ""), 6, 5)
  print(pt)
  dev.off()
  
  # relative performance output data
  fwrite(dt_one_cor, file = paste(result_dir, "/performance_assessment/performance_evaluation_in_relative_expression.txt", sep = ""), sep = "\t")
  
  }