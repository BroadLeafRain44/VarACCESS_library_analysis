#!/usr/bin/env Rscript
#
# Fit log_correct_GAM on neg/pos controls, apply to top50series (no refit).
# Requires clean_outputs/raw_data_read50_min2replicates.tsv.
#
# Model:
#   log(rate + eps) ~ Classification + replicate + s(Total_Editable_Positions)
# Corrected:
#   exp(log(rate + eps) - replicate - s(positions)) - eps
#

library(mgcv)

base_dir <- "/Users/jordanjalbert-ross/Desktop/MAITRISE/LAB/DRYLAB/library/2026_02_09_library/final"
source(file.path(base_dir, "simple_scripts/library_helpers.R"))

eps <- 0.001
out_dir <- file.path(base_dir, "clean_outputs")
fit_path <- file.path(out_dir, "log_correct_GAM_TC_GA_Edit_Rate_fit_controls.rds")
out_controls <- file.path(out_dir, "log_correct_GAM_TC_GA_Edit_Rate_predictions.tsv")
out_top50 <- file.path(out_dir, "log_correct_GAM_TC_GA_Edit_Rate_predictions_top50series.tsv")

keep_cols <- c(
  "Reference_ID", "Classification", "TC_GA_Edit_Rate", "Read_Count",
  "Total_Editable_Positions",
  "Effect_replicate", "Effect_Total_Editable_Positions", "Tech_fitted",
  "Corrected_TC_GA_Edit_Rate", "replicate"
)

num_cols <- c(
  "TC_GA_Edit_Rate", "Read_Count", "Total_Editable_Positions",
  "Effect_replicate", "Effect_Total_Editable_Positions", "Tech_fitted",
  "Corrected_TC_GA_Edit_Rate"
)

correct_rates <- function(dat, fit) {
  terms <- predict(fit, newdata = dat, type = "terms")
  dat$Effect_replicate <- terms[, "replicate"]
  dat$Effect_Total_Editable_Positions <- terms[, "s(Total_Editable_Positions)"]
  dat$Tech_fitted <- dat$Effect_replicate + dat$Effect_Total_Editable_Positions
  dat$Corrected_TC_GA_Edit_Rate <- exp(
    log(dat$TC_GA_Edit_Rate + eps) - dat$Tech_fitted
  ) - eps
  dat$replicate <- as.character(dat$replicate)
  dat
}

write_with_avg <- function(dat, path) {
  dat <- dat[, keep_cols, drop = FALSE]
  avg <- aggregate(
    dat[num_cols],
    by = list(Reference_ID = dat$Reference_ID, Classification = dat$Classification),
    FUN = mean,
    na.rm = TRUE
  )
  avg$replicate <- "avg"
  cat("avg rows:", nrow(avg), "\n")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(rbind(dat, avg[, keep_cols]), path, sep = "\t", row.names = FALSE, quote = FALSE)
}

df_all <- read_all_sheets(filtered_raw_data_path(base_dir))
ok <- stats::complete.cases(df_all[, c("TC_GA_Edit_Rate", "Total_Editable_Positions")])
df_all <- df_all[ok, , drop = FALSE]

df_all$Classification <- classify(df_all$Reference_ID)

# --- fit on controls ---------------------------------------------------------

ctrl <- df_all[df_all$Classification %in% c("negative", "positive"), , drop = FALSE]
ctrl$Classification <- factor(ctrl$Classification, levels = c("negative", "positive"))
ctrl$replicate <- factor(ctrl$replicate, levels = sort(unique(ctrl$replicate)))
cat("Fitting on", nrow(ctrl), "control rows\n")

fit <- gam(
  I(log(TC_GA_Edit_Rate + eps)) ~ Classification + replicate +
    s(Total_Editable_Positions, k = 10),
  data = ctrl,
  method = "REML"
)
saveRDS(fit, fit_path)
print(summary(fit))

write_with_avg(correct_rates(ctrl, fit), out_controls)
cat("Saved:", out_controls, "\n")

# --- apply to top50series ----------------------------------------------------

t50 <- df_all[df_all$Classification == "top50series", , drop = FALSE]
t50$Classification <- factor("positive", levels = c("negative", "positive"))
t50$replicate <- factor(t50$replicate, levels = levels(ctrl$replicate))
cat("Applying fit to", nrow(t50), "top50series rows\n")

t50_out <- correct_rates(t50, fit)
t50_out$Classification <- "top50series"
write_with_avg(t50_out, out_top50)
cat("Saved:", out_top50, "\n")
