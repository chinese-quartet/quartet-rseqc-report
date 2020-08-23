# S1-6 SNR
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

get_pca_list = function(expr_mat_forsignoise, exp_design) {
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

# get lst gene which was used to count snr by reference dataset
get_more_group <- function(exp_fpkm_log, dt_meta, result_dir) {
  exp_design = (dt_meta[, .(library, group = sample)] %>% setkey(., library))
  exp_fpkm_log_fiter <- exp_fpkm_log[ref_data$snr_lst_gene,]
  pca_list <- get_pca_list(exp_fpkm_log_fiter, exp_design)
  fwrite(pca_list, file = paste(result_dir, "/performance_assessment/studydesign_snr.txt", sep = ""), sep = "\t")

  # S3 based on study design
  # S3-1 SD performance summary
  SNR_value <- pca_list$SNR[1]
  dt_test_sd <- data.table(Batch = "Test", SNR = pca_list$SNR[1])
  SD_ref_compare <- rbind(ref_data$SD_ref[, c("Batch", "SNR")], dt_test_sd)

  test_rank <- c(rank(SD_ref_compare$SNR)[which(SD_ref_compare$Batch == "Test")])

  studydesign_performance_summary <- data.table("Term" = "Signal-to-noise ratio (SNR)",
                                                "Value (%)" = SNR_value,
                                                "Reference value" = "[5, inf)",
                                                "Rank" = paste(test_rank, "16", sep = "/")
  )
  fwrite(studydesign_performance_summary, file = paste(result_dir, "/performance_assessment/studydesign_performance_summary_more.txt", sep = ""), sep = "\t")

  # S5 Simplified report output
  # S5-7 SNR figure
  pdf(file = paste(result_dir, "/simplified_report/", "SNR", ".pdf", sep = ""), 6, 5)
  pt <- ggplot(pca_list, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample), size = 2.5) +
    theme_few() +
    make_theme() +
    guides(shape = guide_legend(ncol = 1), color = guide_legend(ncol = 1, title.position = "top")) +
    theme(axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14)) +
    scale_fill_manual(values = c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745")) +
    scale_color_manual(values = c("#2f5c85", "#7BC8A4", "#FFC65D", "#F16745")) +
    ggtitle(paste("SNR: ", pca_list$SNR[1], sep = "")) +
    labs(x = paste("PC1 (", pca_list$PC1_ratio, "%)", sep = ""),
         y = paste("PC1 (", pca_list$PC2_ratio, "%)", sep = ""))
  print(pt)
  dev.off()

  # performance score figure
  # count mean of each batch in JI, CV, CTR and SNR
  SD_performance_mean_more <<- SD_ref_compare
}
