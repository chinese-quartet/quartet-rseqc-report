#' Get one group
#'
#' @importFrom stats cor
#' @importFrom ggplot2 ggtitle
#' @importFrom ggthemes theme_few
#' @importFrom ggplot2 scale_fill_viridis_c
#' @importFrom data.table fwrite
#' @importFrom dplyr %>%
#' @export
#' 
get_one_group <- function(dt_fpkm_log, dt_counts, dt_meta, result_dir) {
  sample_type_list <- dt_meta[['sample']] %>% unique()
  
  dt_abs_cor <- lapply(sample_type_list, function(x){
    if(dim(dt_meta[sample == x])[1] < 2){
      stop(paste(x, 'must contain two replicates', sep = ' '))
    } else {
      replicate_combn <- data.table(combn(dt_meta[sample == x][['group']], 2))
      mat_abs_cor_per <-  lapply(1:dim(replicate_combn)[2], function(y){
        compare_name <- paste(replicate_combn[[y]][1], 'vs', replicate_combn[[y]][2], sep = '')
        sample_id1 <- dt_meta[group == replicate_combn[[y]][1]][['library']]
        sample_id2 <-dt_meta[group == replicate_combn[[y]][2]][['library']]
        gene_list_1 <- dt_counts[dt_counts[[sample_id1]] >= 3][['gene_id']]
        gene_list_2 <- dt_counts[dt_counts[[sample_id2]] >= 3][['gene_id']]
        gene_list_com <- intersect(gene_list_1, gene_list_2)
        cor_vlaue <- cor(dt_fpkm_log[gene_list_com, on = .(gene_id)][[sample_id1]],
                         dt_fpkm_log[gene_list_com, on = .(gene_id)][[sample_id2]])
        return(list(
          compare_name = compare_name,
          abs_cor = cor_vlaue,
          gene_num = length(gene_list_com)
        )) 
      }) %>% rbindlist()
    }
    
    return(
      mat_abs_cor_per
    )
  }) %>% rbindlist()
  
  ### abs expression replicates correlation of median group ----------------------
  dt_abs_cor_median <- dt_abs_cor[with(dt_abs_cor, which.min(abs(dt_abs_cor$abs_cor - median(dt_abs_cor$abs_cor))))]
  output_abs_rep_res <- function(dt_abs_cor, dt_fpkm_log, dt_counts, dt_meta, result_dir){
    compare_name_median <- dt_abs_cor_median[['compare_name']]
    sample_1 <- strsplit(compare_name_median, 'vs')[[1]][1]
    sample_2 <- strsplit(compare_name_median, 'vs')[[1]][2]
    sample_id1 <- dt_meta[group == sample_1][['library']]
    sample_id2 <- dt_meta[group == sample_2][['library']]
    gene_list_1 <- dt_counts[dt_counts[[sample_id1]] >= 3][['gene_id']]
    gene_list_2 <- dt_counts[dt_counts[[sample_id2]] >= 3][['gene_id']]
    gene_list_com <- intersect(gene_list_1, gene_list_2)
    
    dt_one_cor <- dt_fpkm_log[gene_list_com, on = .(gene_id)][, c('gene_id', sample_id1, sample_id2), with = F]
    return(dt_one_cor)
  }
  
  # output abs median cor group scatter plot
  dt_abs_scatter <- output_abs_rep_res(dt_abs_cor, dt_fpkm_log, dt_counts, dt_meta, result_dir)
  colnames(dt_abs_scatter) <- c('gene_id', dt_meta[colnames(dt_abs_scatter)[2:3], on = .(library)][['group']])
  
  # output correlation data
  fwrite(dt_abs_scatter, file = paste(result_dir, "/performance_assessment/absolute_exp_correlation.txt", sep = ""), sep = "\t")
  xlab_abs_cor <- colnames(dt_abs_scatter)[2]
  ylab_abs_cor <- colnames(dt_abs_scatter)[3]
  
  colnames(dt_abs_scatter) <- c('gene_id', 'replicate1', 'replicate2')
  gene_num <- dt_abs_cor_median[['gene_num']]
  cor_vlaue_pt <- round(cor(dt_abs_scatter$replicate1, dt_abs_scatter$replicate2), digits = 3)
  
  # abs correlation plot 
  pt_abs_median_cor <- ggplot2::ggplot(dt_abs_scatter, aes(x = replicate1, y = replicate2)) +
    geom_point(alpha = 0.8, size = 0.3, col = 'steelblue') +
    theme_few() + 
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    scale_fill_viridis_c(name = "density") +
    labs(title = 'One Group: Absolute Correlation',
         subtitle = paste('Correlation = ', cor_vlaue_pt, ' (N = ', gene_num, ')'),
         x = xlab_abs_cor,
         y = ylab_abs_cor)
  
  # return abs meidan correlation value and pt
  # abs_cor_median <- dt_abs_cor_median[['abs_cor']]
  abs_cor_median <- median(dt_abs_cor$abs_cor)
  
  # one group figure and data will be used in more group
  one_group_out <- list(abs_cor_median, pt_abs_median_cor)
  return(one_group_out)
}
