# Purpose:
#   Perform qPCR analysis using the 2^-ΔΔCt method for multiple target genes,
#   using GAPDH as the reference gene and Control as the calibrator group.
#
# Workflow:
#   1. Read raw Ct values
#   2. Reshape data to long format
#   3. Calculate ΔCt, ΔΔCt, and fold change for each target gene
#   4. Perform t-test on ΔCt values between groups
#   5. Save per-gene summary tables
#   6. Generate barplots and boxplots
#
# Input:
#   data/qpcr/qpcr_final_raw_ct_data.csv
#
# Output:
#   results/qpcr/


# Load required package
library(qpcr)
library(dplyr)
library(tidyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

input_file <- "data/qpcr/qpcr_final_raw_ct_data.csv"
output_dir <- "results/qpcr"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ref_gene <- "GAPDH"
calibrator_group <- "Control"
target_genes <- c("ACAN", "CYP2A13", "UGT2B28", "UGT2B11")

get_significance <- function(pval) {
  case_when(
    pval < 0.001 ~ "***",
    pval < 0.01  ~ "**",
    pval < 0.05  ~ "*",
    TRUE         ~ "ns"
  )
}

analyze_gene <- function(df_long, target_gene, ref_gene, calibrator_group, output_dir) {
  message("Analyzing: ", target_gene)

  df_dct <- df_long %>%
    filter(gene %in% c(target_gene, ref_gene)) %>%
    pivot_wider(names_from = gene, values_from = ct) %>%
    filter(!is.na(.data[[target_gene]]), !is.na(.data[[ref_gene]])) %>%
    mutate(delta_ct = .data[[target_gene]] - .data[[ref_gene]])

  calibrator_mean <- df_dct %>%
    filter(Group == calibrator_group) %>%
    summarise(mean_delta_ct = mean(delta_ct, na.rm = TRUE)) %>%
    pull(mean_delta_ct)

  if (length(calibrator_mean) == 0 || is.na(calibrator_mean)) {
    warning("Skipping ", target_gene, ": no calibrator group values found.")
    return(NULL)
  }

  df_ddct <- df_dct %>%
    mutate(
      ddct = delta_ct - calibrator_mean,
      fold_change = 2^(-ddct)
    )

  ttest <- tryCatch(
    t.test(delta_ct ~ Group, data = df_dct),
    error = function(e) NULL
  )

  pval <- if (!is.null(ttest)) ttest$p.value else NA_real_
  significance <- if (!is.na(pval)) get_significance(pval) else NA_character_

  result_summary <- df_ddct %>%
    group_by(Group) %>%
    summarise(
      mean_ddCt = mean(ddct, na.rm = TRUE),
      mean_FC = mean(fold_change, na.rm = TRUE),
      sd_FC = sd(fold_change, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      gene = target_gene,
      p_value = pval,
      significance = significance
    ) %>%
    select(gene, Group, mean_ddCt, mean_FC, sd_FC, p_value, significance)

  write_csv(result_summary, file.path(output_dir, paste0(target_gene, "_summary.csv")))
  write_csv(df_ddct, file.path(output_dir, paste0(target_gene, "_sample_level.csv")))

  plot_data <- df_ddct %>%
    group_by(Group) %>%
    summarise(
      mean_fc = mean(fold_change, na.rm = TRUE),
      sd_fc = sd(fold_change, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(sig = significance)

  y_bar <- max(plot_data$mean_fc + plot_data$sd_fc, na.rm = TRUE)
  y_box <- max(df_ddct$fold_change, na.rm = TRUE)

  p_bar <- ggplot(plot_data, aes(x = Group, y = mean_fc, fill = Group)) +
    geom_col(color = "black", width = 0.6) +
    geom_errorbar(aes(ymin = mean_fc - sd_fc, ymax = mean_fc + sd_fc), width = 0.2) +
    geom_text(aes(label = sig, y = mean_fc + sd_fc + 0.05 * y_bar), size = 6) +
    labs(
      title = paste("Fold Change of", target_gene, "(2^-ΔΔCt Method)"),
      y = "Fold Change (2^-ΔΔCt)",
      x = NULL
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

  ggsave(
    filename = file.path(output_dir, paste0(target_gene, "_barplot.png")),
    plot = p_bar, width = 6, height = 5, dpi = 300
  )

  p_box <- ggplot(df_ddct, aes(x = Group, y = fold_change, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
    stat_compare_means(
      method = "t.test",
      label = "p.signif",
      label.y = y_box + 0.1 * y_box
    ) +
    labs(
      title = paste("Fold Change of", target_gene, "(2^-ΔΔCt Method)"),
      y = "Fold Change (2^-ΔΔCt)",
      x = NULL
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

  ggsave(
    filename = file.path(output_dir, paste0(target_gene, "_boxplot.png")),
    plot = p_box, width = 6, height = 5, dpi = 300
  )

  return(result_summary)
}

df <- read_csv(input_file, show_col_types = FALSE)

required_cols <- c("Sample", "Group", ref_gene)
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df_long <- df %>%
  pivot_longer(
    cols = -c(Sample, Group),
    names_to = "gene",
    values_to = "ct"
  )

all_results <- map_dfr(target_genes, \(g) analyze_gene(df_long, g, ref_gene, calibrator_group, output_dir))

write_csv(all_results, file.path(output_dir, "All_Genes_qPCR_Summary.csv"))
message("Done.")
