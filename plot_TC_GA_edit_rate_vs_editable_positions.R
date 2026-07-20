#!/usr/bin/env Rscript
#
# Editable positions vs raw / corrected TpC edit rate.
#

library(ggplot2)
library(patchwork)

base_dir <- "/PATH/TO/YOUR/PROJECT"  # EDIT ME
source(file.path(base_dir, "git_scripts/library_helpers.R"))  # EDIT if helpers live elsewhere

pred_path <- file.path(base_dir, "clean_outputs/log_correct_GAM_TC_GA_Edit_Rate_predictions.tsv")
plot_path <- file.path(simple_outputs_path(base_dir), "TC_GA_edit_rate_vs_editable_positions_raw_vs_corrected.png")

df <- read.table(pred_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df <- df[as.character(df$replicate) == "avg", , drop = FALSE]
df <- df[classify(df$Reference_ID) %in% c("negative", "positive"), , drop = FALSE]
df <- df[df$Classification %in% c("negative", "positive"), , drop = FALSE]
df$group <- factor(df$Classification, levels = c("negative", "positive"),
                   labels = c("Negative", "Positive"))
df <- df[
  !is.na(df$Total_Editable_Positions) &
    !is.na(df$TC_GA_Edit_Rate) &
    !is.na(df$Corrected_TC_GA_Edit_Rate),
  ,
  drop = FALSE
]

print_cor <- function(dat, rate_col, label) {
  for (g in levels(dat$group)) {
    sub <- dat[dat$group == g, , drop = FALSE]
    ct <- stats::cor.test(sub$Total_Editable_Positions, sub[[rate_col]], method = "pearson")
    fit <- stats::lm(sub[[rate_col]] ~ sub$Total_Editable_Positions)
    cat(label, g, ": r =", signif(unname(ct$estimate), 3),
        ", slope =", signif(coef(fit)[2], 3),
        ", p =", format(ct$p.value, scientific = TRUE), "\n")
  }
}
print_cor(df, "TC_GA_Edit_Rate", "Raw")
print_cor(df, "Corrected_TC_GA_Edit_Rate", "Corrected")

panel <- function(dat, rate_col, ylab, title) {
  ggplot(dat, aes(x = Total_Editable_Positions, y = .data[[rate_col]], color = group)) +
    geom_point(alpha = 0.55, size = 1.3) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
    labs(x = "Editable positions", y = ylab, title = title, color = NULL) +
    theme_classic(base_size = 12)
}

p <- panel(df, "TC_GA_Edit_Rate", "Raw TpC edit rate (%)", "Before correction") +
  panel(df, "Corrected_TC_GA_Edit_Rate", "Corrected TpC edit rate (%)", "After correction") +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
ggsave(plot_path, p, width = 12, height = 5, dpi = 150)
cat("n =", nrow(df), "\n")
cat("Saved:", plot_path, "\n")
