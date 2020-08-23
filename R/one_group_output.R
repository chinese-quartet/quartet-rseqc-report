get_one_group <- function(dt_exp_annot, dt_pairs, result_dir) {
  # S1-2 detected gene number
  dt_detect <- dt_exp_annot[fpkm > 0.1][, .(n_rep = .N), by = .(gene, sample)]
  dt_detect_stats <- dt_detect[, .N, by = .(sample, n_rep)]

  fwrite(dt_detect_stats, file = paste(result_dir, "/quantification_qc/detected_gene_num.txt", sep = ""), sep = "\t")

  # S1-3 detected gene number JI
  dt_jacard_detect_number <- split(dt_pairs, by = c("library_A", "library_B")) %>% lapply(., function(x) {
    library_A = x$library_A
    library_B = x$library_B

    lst_gene_detect_libraryA <- dt_exp_melt[library == library_A][fpkm > 0.1]$gene
    lst_gene_detect_libraryB <- dt_exp_melt[library == library_B][fpkm > 0.1]$gene

    n_intersect <- intersect(lst_gene_detect_libraryA, lst_gene_detect_libraryB) %>% length()
    n_union <- union(lst_gene_detect_libraryA, lst_gene_detect_libraryB) %>% length()

    return(
      list(library_A = library_A,
           library_B = library_B,
           n_intersect = n_intersect,
           n_union = n_union,
           JacardIndex = n_intersect / n_union
      )
    )
  }) %>% rbindlist()
  dt_jacard_detect_number <- cbind(dt_jacard_detect_number, dt_pairs)
  dt_jacard_detect_number_intra <- dt_jacard_detect_number[.("Intra-sample"), on = .(sample_type)]
  dt_jacard_detect_number_mean <- unique(dt_jacard_detect_number_intra$sampleA) %>% lapply(., function(x) {
    JI_mean = mean(dt_jacard_detect_number_intra[.(x), on = .(sampleA)][["JacardIndex"]])

    return(
      list(
        Sample = x,
        JacardIndex = JI_mean
      )
    )
  }) %>% rbindlist()

  # S1-4 CV
  lst_gene <- dt_exp_annot[fpkm > 0.1, .N, by = gene][N == nrow(dt_meta)]$gene
  dt_cv <- dt_exp_annot[gene %in% lst_gene, .(CV = sd(fpkm) / mean(fpkm) * 100, mean = mean(fpkm), n = .N),
                        by = .(sample, gene)]
  fwrite(dt_cv, file = paste(result_dir, "/quantification_qc/sd_one_group_cv.txt", sep = ""), sep = "\t")
  dt_cv_mean <- unique(dt_cv$sample) %>% lapply(., function(x) {
    cv_mean = mean(dt_cv[.(x), on = .(sample)][["CV"]])

    return(
      list(
        sample = x,
        CV = cv_mean)
    )
  }) %>% rbindlist()

  # S1-5 CTR
  dt_LIR2 <- split(dt_pairs, by = c("library_A", "library_B")) %>% lapply(., function(x) {
    library_A <- x$library_A
    library_B <- x$library_B

    lst_gene_detect_libraryA <- dt_exp_annot[library == library_A][fpkm > 0.1]$gene
    lst_gene_detect_libraryB <- dt_exp_annot[library == library_B][fpkm > 0.1]$gene
    lst_gene_for_IR2 <- intersect(lst_gene_detect_libraryA, lst_gene_detect_libraryB)

    exp_fpkm_log_d <- exp_fpkm_log[lst_gene_for_IR2,]
    dt_exp_melt_d5 <- data.table("GENE_ID" = row.names(exp_fpkm_log_d), exp_fpkm_log_d[, grep("D5", colnames(exp_fpkm_log_d))])
    fwrite(dt_exp_melt_d5, file = paste(result_dir, "/performance_assessment/d5_correlation.txt", sep = ""), sep = "\t")

    return(
      list(
        library_A = library_A,
        library_B = library_B,
        n_gene = length(lst_gene_for_IR2),
        Cor = cor(exp_fpkm_log[lst_gene_for_IR2, library_A],
                  exp_fpkm_log[lst_gene_for_IR2, library_B])
      )
    )

  }) %>% rbindlist()
  dt_LIR2_pairs <- cbind(dt_LIR2, dt_pairs)
  dt_LIR2_pairs_intra <- dt_LIR2_pairs[.("Intra-sample"), on = .(sample_type)]
  dt_LIR2_pairs_mean <- unique(dt_LIR2_pairs$sampleA) %>% lapply(., function(x) {
    CTR_mean <- mean(dt_LIR2_pairs_intra[.(x), on = .(sampleA)][["Cor"]])

    return(
      list(
        Sample = x,
        CTR = CTR_mean
      )
    )
  }) %>% rbindlist()

  # S2 depend on reference dataset
  # compare detected gene performance between new dataset and reference dataset
  # S2-1 detected gene
  ref_data$ref_detected_gene_per$Group <- "reference"
  dt_exp_annot[, check_out := sapply(fpkm, function(x) { ifelse(x >= 0.1, "yes", "no") })]

  detected_gene_performance_test <- unique(dt_exp_annot[["group"]]) %>% lapply(., function(x) {

    pos_test <- dt_exp_annot[.(x, "yes"), on = .(group, check_out)][["gene"]]
    neg_test <- dt_exp_annot[.(x, "no"), on = .(group, check_out)][["gene"]]
    sample_type <- unique(dt_exp_annot[.(x, "yes"), on = .(group, check_out)][["sample"]])
    ref_pos <- ref_data$ref_data[.(sample_type, "Pos_tier1"), on = .(Sample, Type)][["Gene"]]
    ref_neg <- ref_data$ref_data[.(sample_type, "Neg_tier1"), on = .(Sample, Type)][["Gene"]]

    TP <- length(intersect(ref_pos, pos_test))
    TN <- length(intersect(ref_neg, neg_test))
    FN <- length(setdiff(ref_pos, pos_test))
    FP <- length(setdiff(ref_neg, neg_test))

    Sensitivity <- TP / (TP + FN)
    Specificity <- TN / (FP + TN)

    return(
      list(Sensitivity_tier1 = Sensitivity,
           Specificity_tier1 = Specificity,
           File = x,
           Sample = sample_type,
           PR = NA,
           Batch = "Test",
           DataQual = NA,
           Platform = NA,
           Group = "Test")
    )
  }) %>% rbindlist()

  detected_gene_performance_compared <- rbind(ref_data$ref_detected_gene_per, detected_gene_performance_test)
  fwrite(detected_gene_performance_compared, file = paste(result_dir, "/quantification_qc/ref_detected_gene_performance_compared.txt", sep = ""), sep = "\t")

  # S2-2 quatification qc table
  detected_gene_summary <- unique(detected_gene_performance_test[["Sample"]]) %>% lapply(., function(x) {
    sensitivity_mean = mean(detected_gene_performance_test[.(x), on = .(Sample)][["Sensitivity_tier1"]])
    specificity_mean = mean(detected_gene_performance_test[.(x), on = .(Sample)][["Specificity_tier1"]])

    return(
      list(
        Sample = x,
        Sensitivity = sensitivity_mean,
        Reference_Sensitivity_min = 0.96,
        Reference_Sensitivity_max = 1,
        Specificity = specificity_mean,
        Reference_Specificity_min = 0.94,
        Reference_Specificity_max = 1
      ))
  }) %>% rbindlist()
  fwrite(detected_gene_summary, file = paste(result_dir, "/quantification_qc/ref_one_group_performance_summary.txt", sep = ""), sep = "\t")


  # S3-1 SD performance summary
  JI_value <- paste(round(mean(dt_jacard_detect_number_intra$JacardIndex), digits = 2), "%+-%",
                    round(sd(dt_jacard_detect_number_intra$JacardIndex), digits = 2), sep = "")
  CV_value <- paste(round(mean(dt_cv$CV), digits = 2), "%+-%",
                    round(sd(dt_cv$CV), digits = 2), sep = "")
  CTR_value <- paste(round(mean(dt_LIR2_pairs_intra$Cor), digits = 2), "%+-%",
                     round(sd(dt_LIR2_pairs_intra$Cor), digits = 2), sep = "")
  dt_test_sd <- data.table(Batch = "Test", Detect.JI = mean(dt_jacard_detect_number_intra$JacardIndex), CV = mean(dt_cv$CV),
                           CTR = mean(dt_LIR2_pairs_intra$Cor))
  SD_ref_compare <- rbind(ref_data$SD_ref[, !"SNR"], dt_test_sd)

  test_rank <- c(rank(SD_ref_compare$Detect.JI)[which(SD_ref_compare$Batch == "Test")],
                 rank(SD_ref_compare$CV)[which(SD_ref_compare$Batch == "Test")],
                 rank(SD_ref_compare$CTR)[which(SD_ref_compare$Batch == "Test")])

  studydesign_performance_summary <- data.table("Term" = c("Detection Jacard Index (JI)",
                                                           "Coefficient of variation (CV)",
                                                           "Correlation of technical relicates (CTR)"),
                                                "Value (%)" = c(JI_value, CV_value, CTR_value),
                                                "Reference value" = c("88.48 %+-% 0.59", "10.24 %+-% 0.63", "97.77 %+-% 0.25"),
                                                "Rank" = paste(test_rank, "16", sep = "/")
  )
  fwrite(studydesign_performance_summary, file = paste(result_dir, "/performance_assessment/studydesign_performance_summary_one.txt", sep = ""), sep = "\t")

  # S5 Simplified report output
  # S5-8 D5_1 and D5_2 correlation figure
  dt_exp_melt_d5 <- fread(paste(result_dir, "/performance_assessment/d5_correlation.txt", sep = ""))
  pdf(file = paste(result_dir, "/simplified_report/", "cor_d5", ".pdf", sep = ""), 6, 6)
  dt_cor <- ggplot(dt_exp_melt_d5, aes(x = FPKM.D5_1, y = FPKM.D5_2)) + geom_point(color = "#2f5c85") +
    theme_few() +
    make_theme() +
    ggtitle(paste("CTR = ", round(dt_LIR2[.("FPKM.D5_1", "FPKM.D5_2"), on = .(library_A, library_B)][["Cor"]], digits = 2), sep = "")) +
    labs(x = "D5_1",
         y = "D5_2")
  print(dt_cor)
  dev.off()

  # S3 based on study design
  sd_one_group_performance_summary <- cbind(
    dt_detect_stats[.(3), on = .(n_rep)][, c("sample", "N")],
    ref_n_min = 20000,
    ref_n_max = 25000,
    dt_jacard_detect_number_mean[, "JacardIndex"],
    ref_ji_min = 0.8,
    ref_ji_max = 1,
    dt_cv_mean[, "CV"],
    ref_cv_min = 5,
    ref_cv_max = Inf,
    dt_LIR2_pairs_mean[, "CTR"],
    ref_ctr_min = 0.95,
    ref_ctr_max = 1)
  fwrite(sd_one_group_performance_summary, file = paste(result_dir, "/quantification_qc/sd_one_group_performance_summary.txt", sep = ""), sep = "\t")


  # performance score figure
  detected_gene_performance_mean <- lapply(ref_data$highqual_batch, function(x) {
    sensitivity_mean <- mean(detected_gene_performance_compared[.(x), on = .(Batch)][["Sensitivity_tier1"]])
    specificity_mean <- mean(detected_gene_performance_compared[.(x), on = .(Batch)][["Specificity_tier1"]])

    return({
      list(
        Batch = x,
        Detected_gene_sensitivity = sensitivity_mean,
        Detected_gene_specificity = specificity_mean
      )
    })
  }) %>% rbindlist()

  # count mean of each batch in JI, CV, CTR
  SD_performance_mean_one <<- SD_ref_compare
  detected_gene_performance_mean <<- detected_gene_performance_mean
}
