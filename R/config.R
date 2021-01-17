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
  ref_data$qcintra_forplot <- readRDS(paste0(ref_data_dir, "/qcintra_forplot.rds"))
  ref_data$corr_ref <- fread(paste0(ref_data_dir, "/TableS2_ReferenceDatasets.csv"), drop = 'V1')
  ref_data$refqc_202011_forplot <- readRDS(paste0(ref_data_dir, "/refqc_202011_forplot.rds"))

  return(ref_data)
}

#' This function is for identify DEG from a given matrix
#'
#' @param expr_mat should a matrix in stead of a data frame, with column for samples and rows for signature(gene/protein/metabolics...)
#'                 For transcriptome, a matrix after log transformation is recommended.
#' @param group    should be a factor whose length is identical to the number of the columns in expr_mat,
#'                 describing the group information of each column in expr_mat
#' @importFrom edgeR DGEList
#' @importFrom edgeR filterByExpr
#' @importFrom edgeR calcNormFactors
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom stats stats
#' @importFrom dplyr %>%
#' @importFrom data.table as.data.table
#' @export
library(edgeR)
library(limma)

DEGanalysis <- function(exprMat, group){
  dge <- DGEList(counts = exprMat)
  design <- model.matrix(~group)
  
  keep <- filterByExpr(dge, design)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  
  v <- voom(dge, design, plot=F)
  fit <- lmFit(v, design)
  
  fit <- eBayes(fit)
  result <- topTable(fit, coef=ncol(design), sort.by = 'logFC', number = Inf)
  result$gene = rownames(result)
  result$groupA =  levels(group)[1]
  result$groupB =  levels(group)[2]
  return(as.data.table(result))
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


#' Make scatter and box figure
#' @importFrom ggplot2 ggplot
#' @importFrom cowplot

plot_scatter_box <- function(dt_sb, var_x, var_y, col_g, xlab, ylab, title_lab){
  pmain <- ggplot(dt_sb, aes_string(x = var_x, y = var_y, color = col_g)) +
    geom_point() +
    scale_color_manual(values = c('#2f5c85', '#7ba1c7', 'red')) +
    theme_few() +
    theme(legend.position = "none") +
    labs(title = title_lab, x = xlab, y = ylab)
  
  xplot <- ggplot(dt_sb, aes_string(x = col_g, y = var_x, colour = col_g)) + 
    geom_boxplot() +
    scale_color_manual(values = c('#2f5c85', '#7ba1c7', 'red')) +
    coord_flip() +
    theme_classic()
  
  yplot <- ggplot(dt_sb, aes_string(x = col_g, y = var_y, colour = col_g)) + 
    geom_boxplot() +
    scale_color_manual(values = c('#2f5c85', '#7ba1c7', 'red')) +
    theme_classic()
  
  p1 <- insert_xaxis_grob(pmain, xplot, grid::unit(.2, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, yplot, grid::unit(.2, "null"), position = "right")
  pt_sb <- ggdraw(p2)
  return(pt_sb)
}

#' S1-6 SNR
#'
#' @importFrom data.table setkey
calc_signoise_ratio <- function(pca_prcomp, exp_design) {
  
  pcs <- as.data.frame(predict(pca_prcomp))
  dt_perc_pcs <- data.table(PCX = 1:nrow(pcs),
                            Percent = summary(pca_prcomp)$importance[2,],
                            AccumPercent = summary(pca_prcomp)$importance[3,])
  
  dt_dist <- data.table(ID.A = rep(rownames(pcs), each = nrow(pcs)),
                        ID.B = rep(rownames(pcs), time = nrow(pcs)))
  
  dt_dist$Group.A <- exp_design[dt_dist$ID.A]$group
  dt_dist$Group.B <- exp_design[dt_dist$ID.B]$group
  
  dt_dist[, Type := ifelse(ID.A == ID.B, "Same",
                           ifelse(Group.A == Group.B, "Intra", "Inter"))]
  dt_dist[, Dist := sqrt(dt_perc_pcs[1]$Percent * (pcs[ID.A, 1] - pcs[ID.B, 1]) ^ 2 + dt_perc_pcs[2]$Percent * (pcs[ID.A, 2] - pcs[ID.B, 2]) ^ 2)]
  
  dt_dist_stats <- dt_dist[, .(Avg.Dist = mean(Dist)), by = .(Type)]
  setkey(dt_dist_stats, Type)
  signoise <- dt_dist_stats["Inter"]$Avg.Dist / dt_dist_stats["Intra"]$Avg.Dist
  return(signoise)
}

#' Get PCA list
#'
#' @importFrom stats prcomp
get_pca_list <- function(expr_mat_forsignoise, exp_design, dt_meta) {
  pca_prcomp = prcomp(t(expr_mat_forsignoise), scale = F)
  pcs = predict(pca_prcomp) %>% data.frame()
  pcs$library = row.names(pcs)
  pcs_add_meta = merge(pcs, dt_meta, by = "library")
  PC1_ratio = round(summary(pca_prcomp)$importance[2, 1] * 100, digits = 2)
  PC2_ratio = round(summary(pca_prcomp)$importance[2, 2] * 100, digits = 2)
  PC3_ratio = round(summary(pca_prcomp)$importance[2, 3] * 100, digits = 2)
  SNR = round(calc_signoise_ratio(pca_prcomp, exp_design = exp_design), digits = 2)
  gene_num = dim(expr_mat_forsignoise)[1]
  pca_list = cbind(pcs_add_meta, PC1_ratio, PC2_ratio, PC3_ratio, SNR, gene_num)
  return(pca_list)
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