#!/usr/bin/env Rscript
#
# Representative within-series scatters: top / median / bottom by raw r,
# raw (row 1) and LFC (row 2).
#

library(ggplot2)
library(patchwork)

base_dir <- "/Users/jordanjalbert-ross/Desktop/MAITRISE/LAB/DRYLAB/library/2026_02_09_library/final"
source(file.path(base_dir, "simple_scripts/library_helpers.R"))

pred_top50_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions_top50series.tsv")
pred_ctrl_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions.tsv")
dnase_path <- file.path(base_dir, "data/top50series_dnase_scores_per_allele_16KB.tsv")
seq_path <- file.path(base_dir, "data/varACCESS_DddA_firstlibrary.xlsx")
plot_path <- file.path(simple_outputs_path(base_dir), "top50series_alphaG_vs_edit_within_series_representative_raw_lfc.png")

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
  out[order(out$series_idx), , drop = FALSE]
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

ranked <- series_stats_raw[order(-series_stats_raw$pearson_r), , drop = FALSE]
sid_top <- ranked$series_idx[1]
sid_bottom <- ranked$series_idx[nrow(ranked)]
median_r <- median(series_stats_raw$pearson_r, na.rm = TRUE)
d_med <- abs(series_stats_raw$pearson_r - median_r)
ord_med <- order(d_med, series_stats_raw$series_idx)
sid_median <- series_stats_raw$series_idx[ord_med[1]]
chosen <- c(sid_top, sid_median, sid_bottom)

r_of <- function(stats_df, sid) stats_df$pearson_r[stats_df$series_idx == sid]

scatter <- function(dat, x_col, y_col, xlab, ylab, title) {
  ggplot(dat, aes(x = .data[[x_col]], y = .data[[y_col]], color = mutations_nt)) +
    geom_point(alpha = 0.85, size = 2) +
    scale_color_viridis_c(option = "plasma", name = "Mutated bases (nt)") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "black") +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")
}

panels <- list()
for (sid in chosen) {
  panels[[length(panels) + 1L]] <- scatter(
    df_raw[df_raw$series_idx == sid, , drop = FALSE],
    "dnase_score", "edit_rate",
    "AlphaGenome DNase score", "TpC edit rate (%)",
    sprintf("Series %d", sid)
  )
}
for (sid in chosen) {
  panels[[length(panels) + 1L]] <- scatter(
    df_lfc[df_lfc$series_idx == sid, , drop = FALSE],
    "dnase_lfc", "edit_lfc",
    "AlphaGenome DNase LFC", "TpC edit rate LFC",
    sprintf("Series %d (LFC)", sid)
  )
}

p <- wrap_plots(panels, ncol = 3, nrow = 2)

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
ggsave(plot_path, p, width = 12, height = 8, dpi = 150)

for (sid in chosen) {
  cat("Series", sid, ": raw r =", signif(r_of(series_stats_raw, sid), 3),
      ", LFC r =", signif(r_of(series_stats_lfc, sid), 3), "\n")
}
cat("Chosen (by raw r): top =", sid_top, "| median =", sid_median, "| bottom =", sid_bottom, "\n")
cat("Saved:", plot_path, "\n")
