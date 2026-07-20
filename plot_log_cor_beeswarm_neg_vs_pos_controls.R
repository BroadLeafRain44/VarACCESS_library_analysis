#!/usr/bin/env Rscript
#
# Negative vs positive controls: log-corrected TpC edit rate.
#

library(ggplot2)

base_dir <- "/PATH/TO/YOUR/PROJECT"  # EDIT ME
source(file.path(base_dir, "git_scripts/library_helpers.R"))  # EDIT if helpers live elsewhere

pred_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions.tsv")
plot_path <- file.path(simple_outputs_path(base_dir), "log_cor_beeswarm_corrected_TC_GA_neg_vs_pos_controls.png")

df <- read.table(pred_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df <- df[as.character(df$replicate) == "avg", , drop = FALSE]
df <- df[df$Classification %in% c("negative", "positive"), , drop = FALSE]
df$group <- factor(df$Classification, levels = c("negative", "positive"),
                   labels = c("Negative", "Positive"))
df$y <- pmax(0, df$Corrected_TC_GA_Edit_Rate)

tt <- stats::t.test(y ~ group, data = df, var.equal = FALSE)
cat("t-test p =", format(tt$p.value, digits = 3, scientific = TRUE), "\n")

p <- ggplot(df, aes(x = group, y = y)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  labs(x = NULL, y = "TpC edit rate (%)") +
  theme_classic(base_size = 12)

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
ggsave(plot_path, p, width = 5, height = 5, dpi = 150)
cat("Saved:", plot_path, "\n")
