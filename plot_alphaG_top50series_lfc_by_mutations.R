#!/usr/bin/env Rscript
#
# AlphaGenome DNase LFC vs cumulative mutations — trend + boxplot.
# 0 nt = 0 LFC; top50series pos 1..18 = raw_score_log2.
#

library(ggplot2)

base_dir <- "/PATH/TO/YOUR/PROJECT"  # EDIT ME
source(file.path(base_dir, "git_scripts/library_helpers.R"))  # EDIT if helpers live elsewhere

top50_dnase_path <- file.path(base_dir, "data/top50series_dnase_scores_per_allele_16KB.tsv")
seq_path <- file.path(base_dir, "data/varACCESS_DddA_firstlibrary.xlsx")
out_dir <- simple_outputs_path(base_dir)
plot_trend_path <- file.path(out_dir, "alphaG_lfc_by_mutations_trend.png")
plot_box_path <- file.path(out_dir, "alphaG_lfc_by_mutations_boxplot.png")
stats_path <- file.path(out_dir, "alphaG_lfc_by_mutations_pairwise_wilcox.tsv")

top50_map <- build_top50_pos_map(seq_path)
series_ok <- series_with_pos_ctrl_anchor(filtered_raw_data_path(base_dir), top50_map)

top50_dnase <- read.delim(top50_dnase_path, stringsAsFactors = FALSE)
meta <- parse_top50_ids(top50_dnase$Reference_ID)
df <- cbind(top50_dnase, meta[, c("series_idx", "pos_in_series", "mutations_nt")])
df <- df[classify(df$Reference_ID) == "top50series" & df$series_idx %in% series_ok, ]

top50_plot <- data.frame(
  series_idx = df$series_idx, mutations_nt = df$mutations_nt, lfc_plot = df$raw_score_log2
)
ref_zero <- data.frame(
  series_idx = unique(top50_plot$series_idx), mutations_nt = 0L, lfc_plot = 0
)
df_plot <- rbind(ref_zero, top50_plot)

pairwise_stats <- compute_pairwise_mutation_wilcox(df_plot, "lfc_plot")
plots <- mutation_trend_box_plots(df_plot, "lfc_plot", "Log2 FC AlphaGenome DNase", hline0 = TRUE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(plot_trend_path, plots$trend, width = 8, height = 5.5, dpi = 150)
ggsave(plot_box_path, plots$box, width = 10, height = 6, dpi = 150)
write.table(pairwise_stats, stats_path, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Series plotted:", length(unique(df_plot$series_idx)), "\n")
cat("Significant adjacent pairs (BH-adjusted p < 0.05):\n")
print(pairwise_stats[pairwise_stats$stars != "", , drop = FALSE])
cat("Saved:", plot_trend_path, "\n")
cat("Saved:", plot_box_path, "\n")
