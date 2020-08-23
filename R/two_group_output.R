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
#' @importFrom ggthemes theme_few
make_performance_plot <- function(dt_performance, xname = xname, yname = yname, gname = gname, file_name = file_name) {
  dt_plot <- dt_performance[, c(xname, yname, gname), with = F]
  setnames(dt_plot, xname, "xlab")
  setnames(dt_plot, yname, "ylab")
  setnames(dt_plot, gname, "glab")

  # plot figure
  pt_deg <- ggplot(dt_plot, aes(x = xlab, y = ylab)) +
    geom_point(aes(color = glab), size = 2.5) +
    theme_few() +
    make_theme() +
    guides(shape = guide_legend(ncol = 1), color = guide_legend(ncol = 1, title.position = "top")) +
    theme(axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14)) +
    scale_fill_manual(values = c("#2f5c85", "#7ba1c7", "red")) +
    scale_color_manual(values = c("#2f5c85", "#7ba1c7", "red")) +
    labs(x = xname, y = yname, colour = gname)

  xdens <- ggplot(dt_plot, aes(x = glab, y = xlab, colour = glab)) +
    geom_boxplot() +
    scale_color_manual(values = c("#2f5c85", "#7ba1c7", "red")) +
    coord_flip()

  ydens <- ggplot(dt_plot, aes(x = glab, y = ylab, colour = glab)) +
    geom_boxplot() +
    scale_color_manual(values = c("#2f5c85", "#7ba1c7", "red")) +
    theme_classic()
  pt_deg_x <- insert_xaxis_grob(pt_deg, xdens, grid::unit(.2, "null"), position = "top")
  pt_deg_y <- insert_yaxis_grob(pt_deg_x, ydens, grid::unit(.2, "null"), position = "right")

  pdf(file = paste(result_dir, "/simplified_report/", file_name, ".pdf", sep = ""), 6, 5)
  pt <- ggdraw(pt_deg_y)
  print(pt)
  dev.off()
}


get_two_group <- function(dt_exp_annot, dt_pairs, result_dir) {
  ### S2-3 degs performance-------------------------


  lst_sample_all <- dt_meta$sample %>% unique()
  pairs_sample <- combn(lst_sample_all, 2) %>% t %>% as.data.table()
  setnames(pairs_sample, 1:2, c("sampleA", "sampleB"))
  pairs_sample$compare_group = paste(pairs_sample$sampleB, pairs_sample$sampleA, sep = "/")

  dt_DEG <- apply(pairs_sample, 1, function(x) {
    sampleA = x["sampleA"]
    sampleB = x["sampleB"]

    lst_library_sampleA <- dt_meta[sample == sampleA]$library
    lst_library_sampleB <- dt_meta[sample == sampleB]$library
    lst_library_forDEG <- c(lst_library_sampleA, lst_library_sampleB)

    x_DEG <- DEG_analysis(exp_fpkm_log[, lst_library_forDEG],
                         factor(dt_meta[lst_library_forDEG, on = "library"]$sample))
    x_DEG$compare_group = sprintf("%s/%s", sampleB, sampleA)
    return(x_DEG)
  }) %>% rbindlist()

  degs_performance_test <- unique(pairs_sample[["compare_group"]]) %>% lapply(., function(x) {
    print(x)
    pos_test = dt_DEG[.(x, c("Up-regulated", "Down-regulated")), on = .(compare_group, DEGtype)][["gene"]]
    neg_test = dt_DEG[.(x, "not DEG"), on = .(compare_group, DEGtype)][["gene"]]
    ref_pos = ref_data$ref_degs[.(x, "deg"), on = .(compare_group, checkout)][["gene"]]
    ref_neg = ref_data$ref_degs[.(x, "non-deg"), on = .(compare_group, checkout)][["gene"]]

    TP = length(intersect(ref_pos, pos_test))
    TN = length(intersect(ref_neg, neg_test))
    FN = length(setdiff(ref_pos, pos_test))
    FP = length(setdiff(ref_neg, neg_test))

    Sensitivity = TP / (TP + FN)
    Specificity = TN / (FP + TN)

    return(
      list(Sensitivity = Sensitivity,
           Specificity = Specificity,
           Batch = "Test",
           Compare_group = x,
           DataQual = NA,
           Group = "Test")
    )
  }) %>% rbindlist()

  ref_degs_per_m <- ref_data$ref_degs_per[unique(ref_data$ref_meta[, c("Batch", "symbol.libPrep")]), on = "Batch"]
  setnames(ref_degs_per_m, "symbol.libPrep", "Group")
  ref_degs_performance_compared <- rbind(ref_degs_per_m[DataQual == "HighQual"], degs_performance_test)
  fwrite(ref_degs_performance_compared, file = paste(result_dir, "/performance_assessment/ref_degs_performance_compared.txt", sep = ""), sep = "\t")

  # deg performance summary
  ref_degs_performance_compared_mean <- unique(ref_degs_performance_compared$Batch) %>% lapply(., function(x) {
    Sensitivity_mean = mean(ref_degs_performance_compared[.(x), on = .(Batch)][["Sensitivity"]])
    Specificity_mean = mean(ref_degs_performance_compared[.(x), on = .(Batch)][["Specificity"]])
    Group = unique(ref_degs_performance_compared[.(x), on = .(Batch)][["Group"]])

    return({
      list(
        Batch = x,
        Sensitivity_mean = Sensitivity_mean,
        Specificity_mean = Specificity_mean,
        Group = Group)
    })
  }) %>% rbindlist()
  all_rank_senstivity = rank(ref_degs_performance_compared_mean$Sensitivity_mean)[which(ref_degs_performance_compared_mean$Group == "Test")]
  all_rank_specificity = rank(ref_degs_performance_compared_mean$Specificity_mean)[which(ref_degs_performance_compared_mean$Group == "Test")]
  R_rank_senstivity = rank(ref_degs_performance_compared_mean[Group != "P"][["Sensitivity_mean"]])[which(ref_degs_performance_compared_mean[Group != "P"][["Group"]] == "Test")]
  R_rank_specificity = rank(ref_degs_performance_compared_mean[Group != "P"][["Specificity_mean"]])[which(ref_degs_performance_compared_mean[Group != "P"][["Group"]] == "Test")]
  P_rank_senstivity = rank(ref_degs_performance_compared_mean[Group != "R"][["Sensitivity_mean"]])[which(ref_degs_performance_compared_mean[Group != "R"][["Group"]] == "Test")]
  P_rank_specificity = rank(ref_degs_performance_compared_mean[Group != "R"][["Specificity_mean"]])[which(ref_degs_performance_compared_mean[Group != "R"][["Group"]] == "Test")]

  deg_performance_summary <- data.table("Term" = c("DEGs", "All Rank", "RiboZ Rank", "PolyA Rank"),
                                        "Sensitivity (%)" = c(paste(round(ref_degs_performance_compared_mean[.("Test"), on = .(Batch)][["Sensitivity_mean"]], digits = 2), "%+-%", round(sd(ref_degs_performance_compared[.("Test"), on = .(Batch)][["Sensitivity"]]), digits = 2), sep = ""),
                                                              paste(all_rank_senstivity, "/16", sep = ""), paste(R_rank_senstivity, "/11", sep = ""), paste(P_rank_senstivity, "/6", sep = "")),
                                        "Specificity (%)" = c(paste(round(ref_degs_performance_compared_mean[.("Test"), on = .(Batch)][["Specificity_mean"]], digits = 2), "%+-%", round(sd(ref_degs_performance_compared[.("Test"), on = .(Batch)][["Specificity"]]), digits = 2), sep = ""),
                                                              paste(all_rank_senstivity, "/16", sep = ""), paste(R_rank_specificity, "/11", sep = ""), paste(P_rank_specificity, "/6", sep = "")))

  fwrite(deg_performance_summary, file = paste(result_dir, "/performance_assessment/deg_performance_summary.txt", sep = ""), sep = "\t")

  ### S2-4 reletive expression--------------------
  rel_exp_per_test <- unique(dt_DEG[["compare_group"]]) %>% lapply(., function(x) {
    int_gene = intersect(dt_DEG[x, on = .(compare_group)][["gene"]], ref_data$ref_range[x, on = .(compare_group)][["gene"]])
    dt_DEG_int = dt_DEG[.(x, int_gene), on = .(compare_group, gene)]
    ref_range_int = ref_data$ref_range[.(x, int_gene), on = .(compare_group, gene)]

    return(
      list(
        batch_pairs = "",
        consistent = length(which(dt_DEG_int$logFC >= ref_range_int$ref_lsd & dt_DEG_int$logFC)) / nrow(dt_DEG_int),
        Batch = "Test",
        pairs.x = x,
        corr = cor(dt_DEG_int$logFC, ref_range_int$logFC),
        batch.y = "Test",
        pairs.y = x,
        DataQual = "",
        seq_pro = "",
        Group = "Test"
      ))
  }) %>% rbindlist()
  setnames(ref_data$ref_rel_exp_per, "batch.x", "Batch")
  ref_degs_per_m <- ref_data$ref_rel_exp_per[unique(ref_data$ref_meta[, c("Batch", "symbol.libPrep")]), on = "Batch"]
  setnames(ref_degs_per_m, "symbol.libPrep", "Group")
  ref_rel_exp_per_compared <- rbind(ref_degs_per_m, rel_exp_per_test)
  fwrite(ref_rel_exp_per_compared, file = paste(result_dir, "/performance_assessment/ref_rel_exp_per_compared.txt", sep = ""), sep = "\t")

  # relative expression performance summary
  ref_rel_exp_per_compared_mean <- unique(ref_rel_exp_per_compared$Batch) %>% lapply(., function(x) {
    consistent_mean = mean(ref_rel_exp_per_compared[.(x), on = .(Batch)][["consistent"]])
    corr_mean = mean(ref_rel_exp_per_compared[.(x), on = .(Batch)][["corr"]])
    Group = unique(ref_rel_exp_per_compared[.(x), on = .(Batch)][["Group"]])

    return({
      list(
        Batch = x,
        consistent_mean = consistent_mean,
        corr_mean = corr_mean,
        Group = Group)
    })
  }) %>% rbindlist()

  all_rank_consistent = rank(ref_rel_exp_per_compared_mean$consistent_mean)[which(ref_rel_exp_per_compared_mean$Group == "Test")]
  all_rank_corr = rank(ref_rel_exp_per_compared_mean$corr_mean)[which(ref_rel_exp_per_compared_mean$Group == "Test")]
  R_rank_consistent = rank(ref_rel_exp_per_compared_mean[Group != "P"][["consistent_mean"]])[which(ref_rel_exp_per_compared_mean[Group != "P"][["Group"]] == "Test")]
  R_rank_corr = rank(ref_rel_exp_per_compared_mean[Group != "P"][["corr_mean"]])[which(ref_rel_exp_per_compared_mean[Group != "P"][["Group"]] == "Test")]
  P_rank_consistent = rank(ref_rel_exp_per_compared_mean[Group != "R"][["consistent_mean"]])[which(ref_rel_exp_per_compared_mean[Group != "R"][["Group"]] == "Test")]
  P_rank_corr = rank(ref_rel_exp_per_compared_mean[Group != "R"][["corr_mean"]])[which(ref_rel_exp_per_compared_mean[Group != "R"][["Group"]] == "Test")]

  rel_performance_summary <- data.table("Term" = c("Relative Expression", "All Rank", "RiboZ Rank", "PolyA Rank"),
                                        "Consistent (%)" = c(paste(round(ref_rel_exp_per_compared_mean[.("Test"), on = .(Batch)][["consistent_mean"]], digits = 2), "%+-%", round(sd(ref_rel_exp_per_compared[.("Test"), on = .(Batch)][["consistent"]]), digits = 2), sep = ""),
                                                             paste(all_rank_consistent, "/16", sep = ""), paste(R_rank_consistent, "/11", sep = ""), paste(P_rank_consistent, "/6", sep = "")),
                                        "Correlation (%)" = c(paste(round(ref_rel_exp_per_compared_mean[.("Test"), on = .(Batch)][["corr_mean"]], digits = 2), "%+-%", round(sd(ref_rel_exp_per_compared[.("Test"), on = .(Batch)][["consistent"]]), digits = 2), sep = ""),
                                                              paste(all_rank_corr, "/16", sep = ""), paste(R_rank_corr, "/11", sep = ""), paste(P_rank_corr, "/6", sep = "")))

  fwrite(rel_performance_summary, file = paste(result_dir, "/performance_assessment/rel_performance_summary.txt", sep = ""), sep = "\t")


  ############################# S5 Simplified report output #################################
  # two sample performance table which contain relative expression and degs
  df_two_group_performance_summary <- cbind(rel_exp_per_test[, c("pairs.x", "consistent", "corr")],
                                            consistent_min = 0.82, consistent_max = 1, corr_min = 0.96, corr_max = 1,
                                            degs_performance_test[, c("Sensitivity", "Specificity")],
                                            sensitivity_min = 0.80, sensitivity_max = 1, specificity_min = 0.95, specificity_max = 1)
  colnames(df_two_group_performance_summary)[1] <- "compared_group"
  fwrite(df_two_group_performance_summary, file = paste(result_dir, "/quantification_qc/ref_two_group_performance_summary.txt", sep = ""), sep = "\t")

  # S5-2 DEG performance summary table 
  fwrite(deg_performance_summary, file = paste(result_dir, "/simplified_report/DEG_performace_table.txt", sep = ""), sep = "\t")

  # S5-3 DEG performance figure
  make_performance_plot(dt_performance = ref_degs_performance_compared, xname = "Specificity", yname = "Sensitivity", gname = "Group", file_name = "DEG_per_fig")

  # S5-4 relative expression summary table
  fwrite(rel_performance_summary, file = paste(result_dir, "/simplified_report/DEG_performace_table.txt", sep = ""), sep = "\t")

  # S5-5 relative expression figure
  make_performance_plot(dt_performance = ref_rel_exp_per_compared, xname = "corr", yname = "consistent", gname = "Group", file_name = "rel_per_fig")

  ### count performance score figure ----------------------------------
  # count mean of each batch in relative expression
  rel_exp_performance_mean <<- ref_rel_exp_per_compared_mean[, !"Group"]

  # count mean of each batch in DEGs
  degs_performance_mean <<- ref_degs_performance_compared_mean[, !"Group"]
}