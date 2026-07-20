#!/usr/bin/env Rscript
# Filter complete_results Excel → clean_outputs/raw_data_read50_min2replicates.tsv
# Keep rows with Read_Count >= 50, and only IDs with >= 2 such replicates (samples 1–3).

base_dir <- "/PATH/TO/YOUR/PROJECT"  # EDIT ME
source(file.path(base_dir, "git_scripts/library_helpers.R"))  # EDIT if helpers live elsewhere

xlsx_path <- file.path(base_dir, "data/complete_results_2026_02_09.xlsx")
out_path <- filtered_raw_data_path(base_dir)
min_reps_per_id <- 2L

sheet_map <- data.frame(
  sheet = c("sample1_main", "sample2_main", "sample3_main"),
  replicate = 1:3,
  stringsAsFactors = FALSE
)

df <- read_all_sheets_full_from_excel(xlsx_path, sheet_map)
n_all <- nrow(df)

df <- df[!is.na(df$Read_Count) & df$Read_Count >= FILTER_MIN_READ_COUNT, , drop = FALSE]
n_after_read <- nrow(df)

keep <- names(which(table(df$Reference_ID) >= min_reps_per_id))
df <- df[df$Reference_ID %in% keep, , drop = FALSE]

df$Classification <- classify(df$Reference_ID)
df$group <- df$Classification
df <- df[order(df$Reference_ID, df$replicate), , drop = FALSE]

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.table(df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Rows in:", n_all, "→ after Read_Count >=", FILTER_MIN_READ_COUNT, ":", n_after_read,
    "→ after >=", min_reps_per_id, "reps/ID:", nrow(df),
    "(", length(unique(df$Reference_ID)), "IDs)\n", sep = " ")
cat("Saved:", out_path, "\n")
