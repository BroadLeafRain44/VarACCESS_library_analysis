#!/usr/bin/env Rscript
#
# Within-series Pearson r boxplots (raw and LFC), horizontal.
#

library(ggplot2)
library(patchwork)

base_dir <- "/Users/jordanjalbert-ross/Desktop/MAITRISE/LAB/DRYLAB/library/2026_02_09_library/final"
source(file.path(base_dir, "simple_scripts/library_helpers.R"))

pred_top50_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions_top50series.tsv")
pred_ctrl_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions.tsv")
dnase_path <- file.path(base_dir, "data/top50series_dnase_scores_per_allele_16KB.tsv")
seq_path <- file.path(base_dir, "data/varACCESS_DddA_firstlibrary.xlsx")
out_dir <- simple_outputs_path(base_dir)
boxplot_path <- file.path(out_dir, "top50series_alphaG_vs_edit_within_series_r_boxplot_combined_horizontal.png")
stats_path <- file.path(out_dir, "top50series_alphaG_vs_edit_within_series_r_boxplot_combined_horizontal_stats.tsv")

edit_pseudocount <- 0.01

edit_lfc <- function(rate_pct, baseline_pct, pc = edit_pseudocount) {
  log2((pmax(rate_pct, 0) + pc) / (pmax(baseline_pct, 0) + pc))
}

cor_stats <- function(x, y) {
  ok <- stats::complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3L) return(list(n = length(x), pearson_r = NA_real_, p_value = NA_real_))
  ct <- stats::cor.test(x, y, method = "pearson")
  list(n = length(x), pearson_r = unname(ct$estimate), p_value = ct$p.value)
}

series_r_table <- function(df, x_col, y_col) {
  parts <- lapply(split(df, df$series_idx), function(d) {
    st <- cor_stats(d[[x_col]], d[[y_col]])
    data.frame(series_idx = d$series_idx[1], n = st$n, pearson_r = st$pearson_r,
               p_value = st$p_value, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  out$p_adj_BH <- stats::p.adjust(out$p_value, method = "BH")
  out[order(out$series_idx), , drop = FALSE]
}

h_box <- function(stats_df, xlab) {
  ggplot(stats_df, aes(x = pearson_r, y = 1)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_boxplot(width = 0.4, outlier.shape = NA) +
    geom_jitter(height = 0.15, alpha = 0.7, size = 1.5) +
    scale_y_continuous(breaks = NULL) +
    labs(x = xlab, y = NULL) +
    theme_classic(base_size = 12)
}

top50_map <- build_top50_pos_map(seq_path)
ctrl_avg <- read.table(pred_ctrl_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ctrl_avg <- ctrl_avg[as.character(ctrl_avg$replicate) == "avg", , drop = FALSE]

baseline <- data.frame(
  series_idx = top50_map$series_idx,
  baseline_edit_pct = vapply(top50_map$pos_ctrl_id, mean_corrected_rate, numeric(1), df = ctrl_avg)
)
baseline <- baseline[!is.na(baseline$baseline_edit_pct), , drop = FALSE]

pred_top50 <- read.table(pred_top50_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pred_top50 <- pred_top50[as.character(pred_top50$replicate) == "avg",
                         c("Reference_ID", "Corrected_TC_GA_Edit_Rate"), drop = FALSE]

dnase <- read.table(dnase_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
meta <- parse_top50_ids(dnase$Reference_ID)
dnase <- cbind(dnase, meta[, c("series_idx", "pos_in_series", "mutations_nt")])
dnase_meta <- dnase[
  grepl("^top50series_", dnase$Reference_ID) & dnase$series_idx %in% baseline$series_idx,
  , drop = FALSE
]

mut_raw <- merge(dnase_meta, pred_top50, by = "Reference_ID")
mut_raw <- data.frame(
  series_idx = mut_raw$series_idx,
  mutations_nt = mut_raw$mutations_nt,
  edit_rate = pmax(0, mut_raw$Corrected_TC_GA_Edit_Rate),
  dnase_score = mut_raw$sum_ALT
)
ref_raw <- dnase_meta[dnase_meta$pos_in_series == 1L, , drop = FALSE]
ref_raw <- data.frame(
  series_idx = ref_raw$series_idx, mutations_nt = 0L, dnase_score = ref_raw$sum_REF
)
ref_raw <- merge(ref_raw, baseline, by = "series_idx")
names(ref_raw)[names(ref_raw) == "baseline_edit_pct"] <- "edit_rate"
df_raw <- rbind(mut_raw, ref_raw[, c("series_idx", "mutations_nt", "edit_rate", "dnase_score")])
df_raw <- df_raw[stats::complete.cases(df_raw$dnase_score, df_raw$edit_rate), , drop = FALSE]
series_stats_raw <- series_r_table(df_raw, "dnase_score", "edit_rate")

edit_meta <- parse_top50_ids(pred_top50$Reference_ID)
edit_join <- cbind(pred_top50, edit_meta[, c("series_idx", "mutations_nt")])
edit_join <- merge(edit_join, baseline, by = "series_idx")
edit_join <- edit_join[!is.na(edit_join$mutations_nt), , drop = FALSE]
edit_lfc_df <- data.frame(
  series_idx = edit_join$series_idx,
  mutations_nt = edit_join$mutations_nt,
  edit_lfc = edit_lfc(edit_join$Corrected_TC_GA_Edit_Rate, edit_join$baseline_edit_pct)
)
edit_zero <- data.frame(series_idx = baseline$series_idx, mutations_nt = 0L, edit_lfc = 0)
dnase_lfc_df <- data.frame(
  series_idx = dnase_meta$series_idx[!is.na(dnase_meta$mutations_nt)],
  mutations_nt = dnase_meta$mutations_nt[!is.na(dnase_meta$mutations_nt)],
  dnase_lfc = dnase_meta$raw_score_log2[!is.na(dnase_meta$mutations_nt)]
)
dnase_zero <- data.frame(
  series_idx = unique(dnase_lfc_df$series_idx), mutations_nt = 0L, dnase_lfc = 0
)
df_lfc <- merge(
  rbind(edit_zero, edit_lfc_df),
  rbind(dnase_zero, dnase_lfc_df),
  by = c("series_idx", "mutations_nt")
)
df_lfc <- df_lfc[stats::complete.cases(df_lfc$edit_lfc, df_lfc$dnase_lfc), , drop = FALSE]
series_stats_lfc <- series_r_table(df_lfc, "dnase_lfc", "edit_lfc")

p <- h_box(series_stats_raw, "r (edit rate vs DNase score)") /
  h_box(series_stats_lfc, "r (edit rate LFC vs DNase LFC)")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(boxplot_path, p, width = 10, height = 4.5, dpi = 150)

stats_out <- rbind(
  cbind(series_stats_raw, metric = "raw"),
  cbind(series_stats_lfc, metric = "lfc")
)
write.table(stats_out, stats_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Raw median r =", signif(median(series_stats_raw$pearson_r), 3), "\n")
cat("LFC median r =", signif(median(series_stats_lfc$pearson_r), 3), "\n")
cat("Saved:", boxplot_path, "\n")
