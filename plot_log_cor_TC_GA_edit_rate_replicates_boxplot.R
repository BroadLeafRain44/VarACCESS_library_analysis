#!/usr/bin/env Rscript
#
# Corrected TpC edit rate across replicates.
#

library(ggplot2)
library(dplyr)

base_dir <- "/PATH/TO/YOUR/PROJECT"  # EDIT ME
source(file.path(base_dir, "git_scripts/library_helpers.R"))  # EDIT if helpers live elsewhere

pred_controls <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions.tsv")
pred_top50 <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions_top50series.tsv")
plot_path <- file.path(simple_outputs_path(base_dir), "log_cor_TC_GA_edit_rate_replicates_boxplot.png")

read_pred_replicates <- function(path) {
  d <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  d <- d[!as.character(d$replicate) %in% c("avg", "NA"), , drop = FALSE]
  d$replicate <- as.integer(d$replicate)
  d
}

df <- bind_rows(read_pred_replicates(pred_controls), read_pred_replicates(pred_top50))
df <- df[df$replicate %in% c(1L, 2L, 3L) & !is.na(df$Corrected_TC_GA_Edit_Rate), , drop = FALSE]
df$sample <- factor(paste0("Replicate ", df$replicate),
                    levels = c("Replicate 1", "Replicate 2", "Replicate 3"))

p <- ggplot(df, aes(x = sample, y = Corrected_TC_GA_Edit_Rate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.25, size = 0.7) +
  labs(x = NULL, y = "Corrected TpC edit rate (%)") +
  theme_classic(base_size = 12)

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
ggsave(plot_path, p, width = 6, height = 5, dpi = 150)
cat("Rows plotted:", nrow(df), "\n")
cat("Saved:", plot_path, "\n")
