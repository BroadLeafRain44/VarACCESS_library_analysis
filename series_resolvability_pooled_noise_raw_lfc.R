#!/usr/bin/env Rscript
#
# Per-series resolvability of editing effect above pooled technical noise.
#
# For each series i and replicate j:
#   raw:  d_ij = TpC_ij^Nmut - TpC_ij^0mut
#   LFC:  d_ij = log2((TpC_ij^Nmut + pc) / (TpC_ij^0mut + pc))
#
# Pooled variance: s^2_pool = sum_i sum_j (d_ij - dbar_i)^2 / (S (n - 1))
# One-sample t-test of dbar_i vs SE = sqrt(s^2_pool / n); BH q < 0.05 = resolved.
#

library(ggplot2)

base_dir <- "/PATH/TO/YOUR/PROJECT"  # EDIT ME
source(file.path(base_dir, "git_scripts/library_helpers.R"))  # EDIT if helpers live elsewhere

pred_ctrl_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions.tsv")
pred_top50_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions_top50series.tsv")
seq_path <- file.path(base_dir, "data/varACCESS_DddA_firstlibrary.xlsx")
out_dir <- simple_outputs_path(base_dir)
long_path <- file.path(out_dir, "series_resolvability_pooled_noise_long_table.tsv")
series_path <- file.path(out_dir, "series_resolvability_pooled_noise_per_series.tsv")
summary_path <- file.path(out_dir, "series_resolvability_pooled_noise_summary.tsv")
plot_path_raw <- file.path(out_dir, "series_resolvability_pooled_noise_boxplot_raw.png")
plot_path_lfc <- file.path(out_dir, "series_resolvability_pooled_noise_boxplot_lfc.png")

z_threshold <- 2
mut_levels <- c(10L, 20L, 30L)
scales <- c("raw", "lfc")
edit_pseudocount <- 0.01
bh_alpha <- 0.05
set.seed(1L)

replicate_values <- function(df, reference_id) {
  rows <- df[as.character(df$Reference_ID) == as.character(reference_id) &
               as.character(df$replicate) != "avg", , drop = FALSE]
  if (nrow(rows) == 0L) return(numeric(0))
  v <- rows$Corrected_TC_GA_Edit_Rate
  names(v) <- as.character(rows$replicate)
  v[!is.na(v)]
}

compute_d <- function(df, scale_name) {
  out <- df[!is.na(df$edit_0mut) & !is.na(df$edit_mut), , drop = FALSE]
  if (scale_name == "raw") {
    out$d <- out$edit_mut - out$edit_0mut
  } else {
    out$d <- log2((pmax(out$edit_mut, 0) + edit_pseudocount) /
                    (pmax(out$edit_0mut, 0) + edit_pseudocount))
  }
  out$d_abs <- abs(out$d)
  out$series <- factor(out$series)
  out
}

pooled_variance <- function(df) {
  dbar <- tapply(df$d, df$series, mean)
  resid <- df$d - dbar[df$series]
  S <- length(dbar)
  n <- length(unique(df$rep))
  list(spool2 = sum(resid^2) / (S * (n - 1)), S = S, n = n, df_resid = S * (n - 1))
}

run_resolvability <- function(df, mut_level, scale_name) {
  df <- compute_d(df, scale_name)
  pv <- pooled_variance(df)
  se_mean <- sqrt(pv$spool2 / pv$n)

  series_stats <- do.call(rbind, lapply(split(df, df$series), function(d) {
    data.frame(
      series = as.integer(as.character(d$series[1])),
      d_bar = mean(d$d),
      mean_abs_d = mean(d$d_abs),
      n_rep = nrow(d),
      stringsAsFactors = FALSE
    )
  }))
  series_stats$mut_level <- mut_level
  series_stats$scale <- scale_name
  series_stats$spool2 <- pv$spool2
  series_stats$se_mean <- se_mean
  series_stats$z <- series_stats$d_bar / se_mean
  series_stats$abs_z <- abs(series_stats$z)
  series_stats$pass_z2 <- series_stats$abs_z > z_threshold
  series_stats$pass_z3 <- series_stats$abs_z > 3
  series_stats$t_stat <- series_stats$z
  series_stats$p_value <- 2 * stats::pt(-abs(series_stats$t_stat), df = pv$df_resid)
  series_stats$p_adj_BH <- stats::p.adjust(series_stats$p_value, method = "BH")
  series_stats$resolved <- series_stats$p_adj_BH < bh_alpha

  summary_row <- data.frame(
    mut_level = mut_level, scale = scale_name,
    S = pv$S, n_replicates = pv$n, df_resid = pv$df_resid,
    mean_diff = mean(series_stats$d_bar),
    mean_abs_diff = mean(series_stats$mean_abs_d),
    spool2 = pv$spool2, s_pool = sqrt(pv$spool2), se_mean = se_mean,
    threshold_2se = z_threshold * se_mean,
    n_pass_z2 = sum(series_stats$pass_z2),
    frac_pass_z2 = mean(series_stats$pass_z2),
    n_pass_z3 = sum(series_stats$pass_z3),
    frac_pass_z3 = mean(series_stats$pass_z3),
    n_resolved_BH = sum(series_stats$resolved),
    frac_resolved_BH = mean(series_stats$resolved),
    stringsAsFactors = FALSE
  )
  list(series = series_stats, summary = summary_row)
}

make_plot <- function(series_df, scale_name, ylab, plot_path) {
  plot_df <- series_df[series_df$scale == scale_name, , drop = FALSE]
  plot_df$mut_level <- factor(
    plot_df$mut_level, levels = mut_levels,
    labels = paste0("0 vs ", mut_levels, " mut")
  )
  plot_df$status <- factor(
    ifelse(plot_df$resolved, "Resolved", "Not resolved"),
    levels = c("Resolved", "Not resolved")
  )
  p <- ggplot(plot_df, aes(x = mut_level, y = d_bar, color = status)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_jitter(width = 0.12, alpha = 0.8, size = 1.8) +
    labs(x = NULL, y = ylab, color = NULL) +
    theme_classic(base_size = 12)
  ggsave(plot_path, p, width = 8, height = 4.5, dpi = 150)
}

ctrl_pred <- read.table(pred_ctrl_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
top50_pred <- read.table(pred_top50_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
top50_map <- build_top50_pos_map(seq_path)

long_parts <- lapply(mut_levels, function(m_nt) {
  pos_mut <- as.integer(m_nt / MUTATION_STEP_NT)
  do.call(rbind, lapply(seq_len(nrow(top50_map)), function(i) {
    row <- top50_map[i, , drop = FALSE]
    top50_id_mut <- sprintf(
      "top50series_%d",
      (row$series_idx - 1L) * GUIDES_PER_SERIES + pos_mut
    )
    ref_v <- replicate_values(ctrl_pred, row$pos_ctrl_id)
    mut_v <- replicate_values(top50_pred, top50_id_mut)
    shared <- intersect(names(ref_v), names(mut_v))
    if (length(shared) == 0L) return(NULL)
    data.frame(
      series = row$series_idx, rep = shared, mut_level = m_nt,
      edit_0mut = unname(ref_v[shared]), edit_mut = unname(mut_v[shared]),
      stringsAsFactors = FALSE
    )
  }))
})
long_df <- do.call(rbind, long_parts)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(long_df, long_path, sep = "\t", row.names = FALSE, quote = FALSE)

res_nested <- unlist(lapply(mut_levels, function(m_nt) {
  df_level <- long_df[long_df$mut_level == m_nt, , drop = FALSE]
  lapply(scales, function(sc) run_resolvability(df_level, m_nt, sc))
}), recursive = FALSE)

series_df <- do.call(rbind, lapply(res_nested, `[[`, "series"))
summary_df <- do.call(rbind, lapply(res_nested, `[[`, "summary"))

write.table(series_df, series_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(summary_df, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)

print(summary_df)
for (m_nt in mut_levels) {
  for (sc in scales) {
    level_df <- series_df[series_df$mut_level == m_nt & series_df$scale == sc, , drop = FALSE]
    resolved <- level_df[level_df$resolved, , drop = FALSE]
    resolved <- resolved[order(resolved$series), , drop = FALSE]
    cat(m_nt, "mut,", sc, ":", nrow(resolved), "/", nrow(level_df),
        "resolved (", paste(resolved$series, collapse = ", "), ")\n")
  }
}

make_plot(series_df, "raw", "TpC edit rate difference (%)", plot_path_raw)
make_plot(series_df, "lfc", "TpC edit rate difference (LFC)", plot_path_lfc)
cat("Saved:", plot_path_raw, "\n")
cat("Saved:", plot_path_lfc, "\n")
