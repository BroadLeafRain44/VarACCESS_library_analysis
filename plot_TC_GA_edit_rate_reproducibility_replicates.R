#!/usr/bin/env Rscript
#
# Reproducibility scatters: rep 1 vs 2, 1 vs 3, 2 vs 3 (raw and corrected).
#

library(ggplot2)
library(patchwork)

base_dir <- "/PATH/TO/YOUR/PROJECT"  # EDIT ME
source(file.path(base_dir, "git_scripts/library_helpers.R"))  # EDIT if helpers live elsewhere

out_dir <- simple_outputs_path(base_dir)
plot_raw_path <- file.path(out_dir, "TC_GA_edit_rate_reproducibility_replicates_raw.png")
plot_cor_path <- file.path(out_dir, "log_cor_TC_GA_edit_rate_reproducibility_replicates_corrected.png")
pred_controls <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions.tsv")
pred_top50 <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions_top50series.tsv")

pairs <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L))

pairwise_merge <- function(dat, rate_col, rep_x, rep_y) {
  dx <- dat[dat$replicate == rep_x, c("Reference_ID", "group", rate_col)]
  dy <- dat[dat$replicate == rep_y, c("Reference_ID", rate_col)]
  names(dx)[3] <- "rate_x"
  names(dy)[2] <- "rate_y"
  m <- merge(dx, dy, by = "Reference_ID")
  m[!is.na(m$rate_x) & !is.na(m$rate_y), , drop = FALSE]
}

make_plot <- function(df, rate_col, prefix) {
  panels <- lapply(pairs, function(pr) {
    m <- pairwise_merge(df, rate_col, pr[1], pr[2])
    ggplot(m, aes(x = rate_x, y = rate_y, color = group)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_point(alpha = 0.55, size = 1.3) +
      coord_equal() +
      labs(
        x = paste0(prefix, ", replicate ", pr[1], " (%)"),
        y = paste0(prefix, ", replicate ", pr[2], " (%)"),
        title = paste0("Replicate ", pr[1], " vs ", pr[2]),
        color = NULL
      ) +
      theme_classic(base_size = 12)
  })
  wrap_plots(panels, ncol = 3, guides = "collect") &
    theme(legend.position = "bottom")
}

read_pred_reps <- function(path) {
  d <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  d <- d[!as.character(d$replicate) %in% c("avg", "NA"), , drop = FALSE]
  d$replicate <- as.integer(d$replicate)
  d
}

df_raw <- read_all_sheets(filtered_raw_data_path(base_dir))
df_raw$group <- factor(classify(df_raw$Reference_ID),
                       levels = c("negative", "positive", "top50series"),
                       labels = c("Negative", "Positive", "top50series"))

df_cor <- rbind(read_pred_reps(pred_controls), read_pred_reps(pred_top50))
df_cor$group <- factor(df_cor$Classification,
                       levels = c("negative", "positive", "top50series"),
                       labels = c("Negative", "Positive", "top50series"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(plot_raw_path, make_plot(df_raw, "TC_GA_Edit_Rate", "Raw TpC edit rate"),
       width = 15, height = 5.5, dpi = 150)
ggsave(plot_cor_path, make_plot(df_cor, "Corrected_TC_GA_Edit_Rate", "Corrected TpC edit rate"),
       width = 15, height = 5.5, dpi = 150)
cat("Saved:", plot_raw_path, "\n")
cat("Saved:", plot_cor_path, "\n")
