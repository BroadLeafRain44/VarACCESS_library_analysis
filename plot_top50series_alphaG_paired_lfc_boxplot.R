#!/usr/bin/env Rscript
#
# Paired AlphaGenome DNase LFC: CAD caQTL + top50 0 vs 10/20/30 mut.
#

library(ggplot2)

base_dir <- "/PATH/TO/YOUR/PROJECT"  # EDIT ME
source(file.path(base_dir, "git_scripts/library_helpers.R"))  # EDIT if helpers live elsewhere

top50_dnase_path <- file.path(base_dir, "data/top50series_dnase_scores_per_allele_16KB.tsv")
cad_dnase_path <- file.path(base_dir, "data/CAD_caQTL_dnase_scores_per_allele_16KB.tsv")
out_dir <- simple_outputs_path(base_dir)
cad_selected_path <- file.path(out_dir, "top50series_alphaG_paired_lfc_cad_most_negative50.tsv")
plot_path <- file.path(out_dir, "top50series_alphaG_paired_lfc_boxplot.png")
stats_path <- file.path(out_dir, "top50series_alphaG_paired_lfc_boxplot_stats.tsv")

cad_n <- 50L

top50_dnase <- read.table(top50_dnase_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cad_dnase <- read.table(cad_dnase_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cad_ord <- order(cad_dnase$raw_score_sum, cad_dnase$Reference_ID)
cad_dnase <- cad_dnase[cad_ord, , drop = FALSE]
cad <- data.frame(
  comparison = "CAD",
  lfc = cad_dnase$raw_score_log2[seq_len(cad_n)],
  raw_score_sum = cad_dnase$raw_score_sum[seq_len(cad_n)],
  Reference_ID = cad_dnase$Reference_ID[seq_len(cad_n)],
  stringsAsFactors = FALSE
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(cad, cad_selected_path, sep = "\t", row.names = FALSE, quote = FALSE)

meta <- parse_top50_ids(top50_dnase$Reference_ID)
top50 <- cbind(top50_dnase, meta[, c("series_idx", "pos_in_series")])

paired <- rbind(
  data.frame(comparison = cad$comparison, lfc = cad$lfc, stringsAsFactors = FALSE),
  data.frame(comparison = "0 vs 10", lfc = top50$raw_score_log2[top50$pos_in_series == 1L],
             stringsAsFactors = FALSE),
  data.frame(comparison = "0 vs 20", lfc = top50$raw_score_log2[top50$pos_in_series == 2L],
             stringsAsFactors = FALSE),
  data.frame(comparison = "0 vs 30", lfc = top50$raw_score_log2[top50$pos_in_series == 3L],
             stringsAsFactors = FALSE)
)
paired <- paired[!is.na(paired$lfc), , drop = FALSE]
paired$comparison <- factor(paired$comparison, levels = c("CAD", "0 vs 10", "0 vs 20", "0 vs 30"))

stats_out <- do.call(rbind, lapply(levels(paired$comparison), function(cmp) {
  vals <- paired$lfc[paired$comparison == cmp]
  data.frame(
    comparison = cmp, n = length(vals), mean_lfc = mean(vals),
    p_value = t.test(vals, mu = 0)$p.value, stringsAsFactors = FALSE
  )
}))
write.table(stats_out, stats_path, sep = "\t", row.names = FALSE, quote = FALSE)
print(stats_out)

p <- ggplot(paired, aes(x = comparison, y = lfc)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.6, size = 1.5) +
  labs(x = NULL, y = "AlphaGenome DNase log2 fold-change") +
  theme_classic(base_size = 12)

ggsave(plot_path, p, width = 8, height = 5, dpi = 150)
cat("Saved:", plot_path, "\n")
