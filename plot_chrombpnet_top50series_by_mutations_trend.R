#!/usr/bin/env Rscript
#
# ChromBPNet scores vs cumulative mutations.
#

library(ggplot2)
library(dplyr)
library(readr)

base_dir <- "/Users/jordanjalbert-ross/Desktop/MAITRISE/LAB/DRYLAB/library/2026_02_09_library/final"
source(file.path(base_dir, "simple_scripts/library_helpers.R"))

input_csv <- file.path(base_dir, "data/top_50_series.csv")
output_plot <- file.path(simple_outputs_path(base_dir), "chrombpnet_score_by_mutations_trend.png")

mutation_step_nt <- 10L
max_mutations_nt <- 180L

df <- readr::read_csv(input_csv, show_col_types = FALSE)

loci <- df %>%
  distinct(.data$Chromosome, .data$HG38.start, .data$HG38.end) %>%
  mutate(series_idx = row_number())

df_mut <- df %>%
  left_join(loci, by = c("Chromosome", "HG38.start", "HG38.end")) %>%
  mutate(
    mutations_nt = (as.integer(.data$mutation_round) + 1L) * mutation_step_nt,
    score = .data$Mutant_Score_Avg
  )

wt_rows <- df %>%
  group_by(.data$Chromosome, .data$HG38.start, .data$HG38.end) %>%
  slice(1L) %>%
  ungroup() %>%
  left_join(loci, by = c("Chromosome", "HG38.start", "HG38.end")) %>%
  mutate(mutations_nt = 0L, score = .data$WT_Score_Avg)

df_plot <- bind_rows(wt_rows, df_mut)

trend <- df_plot %>%
  group_by(.data$mutations_nt) %>%
  summarize(y_mean = mean(.data$score, na.rm = TRUE), .groups = "drop")

p <- ggplot(df_plot, aes(x = mutations_nt, y = score)) +
  geom_line(aes(group = series_idx, color = series_idx), alpha = 0.5, linewidth = 0.45) +
  scale_color_viridis_c(option = "plasma", guide = "none") +
  geom_line(data = trend, aes(x = mutations_nt, y = y_mean),
            inherit.aes = FALSE, color = "black", linewidth = 1) +
  labs(x = "Mutated bases (nt)", y = "ChromBPNet score") +
  theme_classic(base_size = 12)

dir.create(dirname(output_plot), recursive = TRUE, showWarnings = FALSE)
ggsave(output_plot, p, width = 7, height = 5, dpi = 150)
cat("Saved:", output_plot, "\n")
