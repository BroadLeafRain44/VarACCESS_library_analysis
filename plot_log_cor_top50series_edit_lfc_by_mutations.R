#!/usr/bin/env Rscript
#
# TpC edit log2 fold change vs cumulative mutations — trend + boxplot.
# LFC = log2((edit_rate + pc) / (baseline + pc)); series with matched pos ctrl.
#

library(ggplot2)

base_dir <- "/PATH/TO/YOUR/PROJECT"  # EDIT ME
source(file.path(base_dir, "git_scripts/library_helpers.R"))  # EDIT if helpers live elsewhere

pred_top50_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions_top50series.tsv")
pred_ctrl_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions.tsv")
seq_path <- file.path(base_dir, "data/varACCESS_DddA_firstlibrary.xlsx")
out_dir <- simple_outputs_path(base_dir)

plot_trend_path <- file.path(out_dir, "log_cor_top50series_edit_lfc_by_mutations_trend.png")
plot_box_path <- file.path(out_dir, "log_cor_top50series_edit_lfc_by_mutations_boxplot.png")
stats_path <- file.path(out_dir, "log_cor_top50series_edit_lfc_by_mutations_pairwise_wilcox.tsv")

edit_pseudocount <- 0.01

edit_lfc <- function(rate_pct, baseline_pct, pseudocount = edit_pseudocount) {
  log2((pmax(rate_pct, 0) + pseudocount) / (pmax(baseline_pct, 0) + pseudocount))
}

top50_map <- build_top50_pos_map(seq_path)
ctrl_avg <- read.table(pred_ctrl_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ctrl_avg <- ctrl_avg[as.character(ctrl_avg$replicate) == "avg", , drop = FALSE]
top50_avg <- read.table(pred_top50_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
top50_avg <- top50_avg[as.character(top50_avg$replicate) == "avg", , drop = FALSE]

baseline <- data.frame(
  series_idx = top50_map$series_idx,
  baseline_edit_pct = vapply(top50_map$pos_ctrl_id, mean_corrected_rate, numeric(1), df = ctrl_avg)
)
baseline <- baseline[!is.na(baseline$baseline_edit_pct), , drop = FALSE]
series_ok <- baseline$series_idx

meta <- parse_top50_ids(top50_avg$Reference_ID)
top50_avg <- cbind(top50_avg, meta[, c("series_idx", "mutations_nt")])
top50_avg <- merge(top50_avg, baseline, by = "series_idx")
top50_avg <- top50_avg[!is.na(top50_avg$mutations_nt) & top50_avg$series_idx %in% series_ok, ]
top50_plot <- data.frame(
  series_idx = top50_avg$series_idx,
  mutations_nt = top50_avg$mutations_nt,
  lfc_plot = edit_lfc(top50_avg$Corrected_TC_GA_Edit_Rate, top50_avg$baseline_edit_pct)
)

ref_zero <- data.frame(
  series_idx = series_ok, mutations_nt = 0L, lfc_plot = 0
)
df_plot <- rbind(ref_zero, top50_plot)

pairwise_stats <- compute_pairwise_mutation_wilcox(df_plot, "lfc_plot")
plots <- mutation_trend_box_plots(df_plot, "lfc_plot", "Log2 FC edit rate", hline0 = TRUE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(plot_trend_path, plots$trend, width = 8, height = 5.5, dpi = 150)
ggsave(plot_box_path, plots$box, width = 10, height = 6, dpi = 150)
write.table(pairwise_stats, stats_path, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Series plotted:", length(unique(df_plot$series_idx)), "\n")
cat("Significant adjacent pairs (BH-adjusted p < 0.05):\n")
print(pairwise_stats[pairwise_stats$stars != "", , drop = FALSE])
cat("Saved trend:", plot_trend_path, "\n")
cat("Saved boxplot:", plot_box_path, "\n")
cat("Saved stats:", stats_path, "\n")
