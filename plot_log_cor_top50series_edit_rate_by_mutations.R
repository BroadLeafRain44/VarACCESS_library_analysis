#!/usr/bin/env Rscript
#
# Log-corrected TpC edit rate vs cumulative mutations — trend + boxplot.
# 0 nt = matched positive control; top50series pos 1..18 = 10..180 nt.
#

library(ggplot2)

base_dir <- "/Users/jordanjalbert-ross/Desktop/MAITRISE/LAB/DRYLAB/library/2026_02_09_library/final"
source(file.path(base_dir, "simple_scripts/library_helpers.R"))

pred_top50_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions_top50series.tsv")
pred_ctrl_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions.tsv")
seq_path <- file.path(base_dir, "data/varACCESS_DddA_firstlibrary.xlsx")
out_dir <- simple_outputs_path(base_dir)
plot_trend_path <- file.path(out_dir, "log_cor_top50series_edit_rate_by_mutations_trend.png")
plot_box_path <- file.path(out_dir, "log_cor_top50series_edit_rate_by_mutations_boxplot.png")
stats_path <- file.path(out_dir, "log_cor_top50series_edit_rate_by_mutations_pairwise_wilcox.tsv")

top50_map <- build_top50_pos_map(seq_path)
ctrl_avg <- read.table(pred_ctrl_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ctrl_avg <- ctrl_avg[as.character(ctrl_avg$replicate) == "avg", , drop = FALSE]
top50_avg <- read.table(pred_top50_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
top50_avg <- top50_avg[as.character(top50_avg$replicate) == "avg", , drop = FALSE]

pos_ctrl_plot <- data.frame(
  series_idx = top50_map$series_idx,
  mutations_nt = 0L,
  edit_rate = vapply(top50_map$pos_ctrl_id, mean_corrected_rate, numeric(1), df = ctrl_avg)
)
pos_ctrl_plot <- pos_ctrl_plot[!is.na(pos_ctrl_plot$edit_rate), , drop = FALSE]
pos_ctrl_plot$edit_rate <- pmax(0, pos_ctrl_plot$edit_rate)
series_ok <- unique(pos_ctrl_plot$series_idx)

meta <- parse_top50_ids(top50_avg$Reference_ID)
top50_avg <- cbind(top50_avg, meta[, c("series_idx", "mutations_nt")])
top50_plot <- top50_avg[
  !is.na(top50_avg$mutations_nt) & top50_avg$series_idx %in% series_ok,
  c("series_idx", "mutations_nt", "Corrected_TC_GA_Edit_Rate")
]
names(top50_plot)[3] <- "edit_rate"
top50_plot$edit_rate <- pmax(0, top50_plot$edit_rate)

df_plot <- rbind(pos_ctrl_plot, top50_plot)
pairwise_stats <- compute_pairwise_mutation_wilcox(df_plot, "edit_rate")
plots <- mutation_trend_box_plots(df_plot, "edit_rate", "TpC edit rate")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(plot_trend_path, plots$trend, width = 8, height = 5.5, dpi = 150)
ggsave(plot_box_path, plots$box, width = 10, height = 6, dpi = 150)
write.table(pairwise_stats, stats_path, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Series plotted:", length(unique(df_plot$series_idx)), "\n")
cat("Significant adjacent pairs (BH-adjusted p < 0.05):\n")
print(pairwise_stats[pairwise_stats$stars != "", , drop = FALSE])
cat("Saved:", plot_trend_path, "\n")
cat("Saved:", plot_box_path, "\n")
