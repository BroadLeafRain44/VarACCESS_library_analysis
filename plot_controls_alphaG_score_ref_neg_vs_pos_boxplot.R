#!/usr/bin/env Rscript
#
# DNase ref scores for negative vs positive controls.
#

library(ggplot2)

base_dir <- "/PATH/TO/YOUR/PROJECT"  # EDIT ME
source(file.path(base_dir, "git_scripts/library_helpers.R"))  # EDIT if helpers live elsewhere

dnase_path <- file.path(base_dir, "data/library_controls_dnase_scores_ref_only.tsv")
filtered_path <- filtered_raw_data_path(base_dir)
plot_path <- file.path(simple_outputs_path(base_dir), "controls_alphaG_score_ref_neg_vs_pos_boxplot.png")

dnase <- read.table(dnase_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
filt <- read.table(filtered_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
library_ids <- unique(filt$Reference_ID)
library_ids <- library_ids[classify(library_ids) %in% c("negative", "positive")]

df <- dnase[dnase$Reference_ID %in% library_ids, , drop = FALSE]
df$Classification <- classify(df$Reference_ID)
df <- df[df$Classification %in% c("negative", "positive"), , drop = FALSE]
df$group <- factor(df$Classification, levels = c("negative", "positive"),
                   labels = c("Negative", "Positive"))
df$y <- df$sum_REF

tt <- stats::t.test(y ~ group, data = df, var.equal = FALSE)
cat("t-test p =", format(tt$p.value, digits = 3, scientific = TRUE), "\n")

p <- ggplot(df, aes(x = group, y = y)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  labs(x = NULL, y = "AlphaGenome DNase score") +
  theme_classic(base_size = 12)

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
ggsave(plot_path, p, width = 5, height = 5, dpi = 150)
cat("Guides plotted:", nrow(df), "\n")
cat("Saved:", plot_path, "\n")
