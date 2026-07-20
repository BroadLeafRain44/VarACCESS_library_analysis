#!/usr/bin/env Rscript
#
# Read count vs GC content with quadratic fit.
#

library(readxl)
library(ggplot2)

base_dir <- "/Users/jordanjalbert-ross/Desktop/MAITRISE/LAB/DRYLAB/library/2026_02_09_library/final"
source(file.path(base_dir, "simple_scripts/library_helpers.R"))

xlsx_path <- file.path(base_dir, "data/complete_results_2026_02_09.xlsx")
seq_path <- file.path(base_dir, "data/varACCESS_DddA_firstlibrary.xlsx")
plot_path <- file.path(simple_outputs_path(base_dir), "read_count_vs_GC_content.png")

gc_pct <- function(seq) {
  s <- toupper(gsub("[^ACGT]", "", as.character(seq)))
  if (!nchar(s)) return(NA_real_)
  100 * sum(strsplit(s, "")[[1]] %in% c("G", "C")) / nchar(s)
}

sheet_map <- data.frame(
  sheet = c("sample1_main", "sample2_main", "sample3_main"),
  replicate = 1:3,
  stringsAsFactors = FALSE
)

va <- as.data.frame(read_excel(seq_path), stringsAsFactors = FALSE)
va$Classification <- classify(va$Sequence_name)
va <- va[va$Classification %in% c("negative", "positive"), , drop = FALSE]
va$seq_190nt <- paste0(va$`9nt_barcode_A1`, va$`181nt_sequence_A1`)
va$GC_pct <- vapply(va$seq_190nt, gc_pct, numeric(1))
va <- va[!is.na(va$GC_pct), c("Sequence_name", "GC_pct", "Classification"), drop = FALSE]
names(va)[1] <- "Reference_ID"

obs_pieces <- vector("list", nrow(sheet_map))
for (i in seq_len(nrow(sheet_map))) {
  d <- as.data.frame(read_excel(xlsx_path, sheet = sheet_map$sheet[i]), stringsAsFactors = FALSE)
  d <- d[, c("Reference_ID", "Read_Count"), drop = FALSE]
  d$replicate <- sheet_map$replicate[i]
  obs_pieces[[i]] <- d
}
obs <- do.call(rbind, obs_pieces)
obs$Classification <- classify(obs$Reference_ID)
obs <- obs[obs$Classification %in% c("negative", "positive"),
  c("Reference_ID", "Read_Count", "replicate", "Classification"), drop = FALSE]

grid <- merge(va, sheet_map[, "replicate", drop = FALSE], by = NULL)
df <- merge(grid, obs, by = c("Reference_ID", "replicate", "Classification"), all.x = TRUE)
df$Read_Count[is.na(df$Read_Count)] <- 0L
df$group <- factor(df$Classification, levels = c("negative", "positive"))
df$log_read <- log10(df$Read_Count + 1)

fit <- lm(log_read ~ GC_pct + I(GC_pct^2), data = df)
r2 <- summary(fit)$r.squared
cat("R2 quadratic =", round(r2, 3), "\n")

xg <- seq(min(df$GC_pct), max(df$GC_pct), length.out = 200)
curve_df <- data.frame(GC_pct = xg, log_read = predict(fit, newdata = data.frame(GC_pct = xg)))
curve_df <- curve_df[curve_df$log_read >= 0, , drop = FALSE]

p <- ggplot(df, aes(GC_pct, log_read, color = group)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_line(data = curve_df, aes(GC_pct, log_read), inherit.aes = FALSE, color = "black", linewidth = 1) +
  scale_color_discrete(labels = c(negative = "Negative", positive = "Positive")) +
  labs(x = "GC (%)", y = "log10(read count + 1)", color = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
ggsave(plot_path, p, width = 6, height = 5, dpi = 150)
cat("Saved:", plot_path, "\n")
