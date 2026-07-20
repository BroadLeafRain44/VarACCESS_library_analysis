#!/usr/bin/env Rscript
#
# DNase score vs TC edit rate (negative / positive / all).
#

library(ggplot2)
library(patchwork)

base_dir <- "/PATH/TO/YOUR/PROJECT"  # EDIT ME
source(file.path(base_dir, "git_scripts/library_helpers.R"))  # EDIT if helpers live elsewhere

dnase_path <- file.path(base_dir, "data/library_controls_dnase_scores_ref_only.tsv")
pred_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions.tsv")
plot_path <- file.path(simple_outputs_path(base_dir), "controls_alphaG_score_vs_TC_edit_rate_correlation.png")

pred <- read.table(pred_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pred <- pred[as.character(pred$replicate) == "avg", , drop = FALSE]
pred <- pred[classify(pred$Reference_ID) %in% c("negative", "positive"), , drop = FALSE]
pred <- pred[pred$Classification %in% c("negative", "positive"), , drop = FALSE]

dnase <- read.table(dnase_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df <- merge(
  pred[, c("Reference_ID", "Classification", "Corrected_TC_GA_Edit_Rate"), drop = FALSE],
  dnase[, c("Reference_ID", "sum_REF"), drop = FALSE],
  by = "Reference_ID"
)
df$group <- factor(df$Classification, levels = c("negative", "positive"),
                   labels = c("Negative", "Positive"))
df <- df[stats::complete.cases(df[, c("sum_REF", "Corrected_TC_GA_Edit_Rate")]), , drop = FALSE]

print_cor <- function(dat, label) {
  ct <- stats::cor.test(dat$sum_REF, dat$Corrected_TC_GA_Edit_Rate, method = "pearson")
  cat(label, ": r =", signif(unname(ct$estimate), 3),
      ", p =", format(ct$p.value, scientific = TRUE), ", n =", nrow(dat), "\n")
}
print_cor(df[df$group == "Negative", , drop = FALSE], "Negative")
print_cor(df[df$group == "Positive", , drop = FALSE], "Positive")
print_cor(df, "All controls")

panel <- function(dat, title, color_by = FALSE) {
  p <- ggplot(dat, aes(x = sum_REF, y = Corrected_TC_GA_Edit_Rate))
  if (color_by) {
    p <- p + geom_point(aes(color = group), alpha = 0.6, size = 1.4)
  } else {
    p <- p + geom_point(alpha = 0.6, size = 1.4)
  }
  p +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, color = "black") +
    labs(x = "AlphaGenome DNase score", y = "TpC edit rate (%)", title = title, color = NULL) +
    theme_classic(base_size = 12)
}

p <- panel(df[df$group == "Negative", , drop = FALSE], "Negative") +
  panel(df[df$group == "Positive", , drop = FALSE], "Positive") +
  panel(df, "All controls", color_by = TRUE) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
ggsave(plot_path, p, width = 14, height = 5, dpi = 150)
cat("Saved:", plot_path, "\n")
