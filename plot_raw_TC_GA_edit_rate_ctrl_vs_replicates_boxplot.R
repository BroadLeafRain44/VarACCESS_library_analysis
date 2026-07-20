#!/usr/bin/env Rscript
#
# Raw TpC edit rate: control vs replicates.
#

library(ggplot2)
library(dplyr)

base_dir <- "/PATH/TO/YOUR/PROJECT"  # EDIT ME
source(file.path(base_dir, "git_scripts/library_helpers.R"))  # EDIT if helpers live elsewhere

xlsx_path <- file.path(base_dir, "data/complete_results_2026_02_09.xlsx")
plot_path <- file.path(simple_outputs_path(base_dir), "raw_TC_GA_edit_rate_ctrl_vs_replicates_boxplot.png")

sheet_map <- data.frame(
  sheet = c("ctrl1_main", "sample1_main", "sample2_main", "sample3_main"),
  replicate = c(0L, 1L, 2L, 3L),
  stringsAsFactors = FALSE
)

df <- read_all_sheets_full_from_excel(xlsx_path, sheet_map)
df$sample <- dplyr::recode(
  df$sheet,
  ctrl1_main = "Control",
  sample1_main = "Replicate 1",
  sample2_main = "Replicate 2",
  sample3_main = "Replicate 3"
)
df$sample <- factor(df$sample, levels = c("Control", "Replicate 1", "Replicate 2", "Replicate 3"))
df <- df %>%
  filter(!is.na(.data$Read_Count), .data$Read_Count >= FILTER_MIN_READ_COUNT, !is.na(.data$TC_GA_Edit_Rate))

p <- ggplot(df, aes(x = sample, y = TC_GA_Edit_Rate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.25, size = 0.7) +
  labs(x = NULL, y = "Raw TpC edit rate (%)") +
  theme_classic(base_size = 12)

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
ggsave(plot_path, p, width = 6, height = 5, dpi = 150)
cat("Rows plotted:", nrow(df), "\n")
cat("Saved:", plot_path, "\n")
