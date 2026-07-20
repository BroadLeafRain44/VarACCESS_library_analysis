#!/usr/bin/env Rscript
#
# ROC/AUC for pos vs neg controls (replicate == "avg").
#

library(ggplot2)

base_dir <- "/Users/jordanjalbert-ross/Desktop/MAITRISE/LAB/DRYLAB/library/2026_02_09_library/final"
source(file.path(base_dir, "simple_scripts/library_helpers.R"))

pred_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions.tsv")
plot_path <- file.path(simple_outputs_path(base_dir), "log_cor_pos_neg_roc_avg_simple.png")

dat <- subset(
  read.delim(pred_path),
  replicate == "avg" & Classification %in% c("negative", "positive") & !is.na(Corrected_TC_GA_Edit_Rate)
)

n0 <- sum(dat$Classification == "negative")
n1 <- sum(dat$Classification == "positive")
wt <- wilcox.test(
  Corrected_TC_GA_Edit_Rate ~ factor(Classification, levels = c("positive", "negative")),
  data = dat, exact = FALSE
)
auc <- as.numeric(wt$statistic) / (n0 * n1)
cat("AUC:", signif(auc, 4), "(n_neg =", n0, ", n_pos =", n1, ")\n")

pos <- as.integer(dat$Classification == "positive")
o <- order(-dat$Corrected_TC_GA_Edit_Rate, -pos)
pos <- pos[o]
curve <- data.frame(
  fpr = c(0, cumsum(pos == 0L) / n0, 1),
  tpr = c(0, cumsum(pos == 1L) / n1, 1)
)

p <- ggplot(curve, aes(fpr, tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_step(direction = "vh") +
  labs(x = "False positive rate", y = "True positive rate",
       title = paste0("AUC = ", signif(auc, 4))) +
  theme_classic(base_size = 12) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
ggsave(plot_path, p, width = 5, height = 5, dpi = 150)
cat("Saved:", plot_path, "\n")
