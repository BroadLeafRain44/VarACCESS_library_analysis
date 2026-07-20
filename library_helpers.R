# Shared data helpers (IDs, paths, Excel/TSV readers, top50series).

FILTERED_RAW_DATA_FILENAME <- "raw_data_read50_min2replicates.tsv"
FILTER_MIN_READ_COUNT <- 50L

GUIDES_PER_SERIES <- 18L
MUTATION_STEP_NT <- 10L
MAX_MUTATIONS_NT <- 180L

classify <- function(reference_id) {
  rid <- as.character(reference_id)
  ifelse(
    grepl("^top50series", rid, ignore.case = TRUE),
    "top50series",
    ifelse(grepl("^Safe", rid, ignore.case = TRUE), "negative", "positive")
  )
}

filtered_raw_data_path <- function(base_dir) {
  file.path(base_dir, paste0("clean_outputs/", FILTERED_RAW_DATA_FILENAME))
}

simple_outputs_path <- function(base_dir) {
  file.path(base_dir, "simple_outputs")
}

read_all_sheets_full_from_excel <- function(path, sheet_map) {
  pieces <- vector("list", nrow(sheet_map))
  for (i in seq_len(nrow(sheet_map))) {
    sh <- sheet_map$sheet[i]
    rep <- sheet_map$replicate[i]
    d <- as.data.frame(readxl::read_excel(path, sheet = sh), stringsAsFactors = FALSE)
    d$replicate <- rep
    d$sheet <- sh
    pieces[[i]] <- d
  }
  do.call(rbind, pieces)
}

# Reads pre-filtered TSV (Read_Count >= 50, >= 2 replicates per ID).
read_all_sheets <- function(path) {
  d <- read.delim(path, stringsAsFactors = FALSE)
  need <- c("Reference_ID", "TC_GA_Edit_Rate", "Read_Count", "Total_Editable_Positions", "replicate")
  d <- d[, need, drop = FALSE]
  d$replicate <- as.integer(d$replicate)
  d
}

# --- top50series helpers -----------------------------------------------------

parse_top50_ids <- function(reference_ids) {
  rid <- as.character(reference_ids)
  n <- as.integer(sub("^top50series_", "", rid, ignore.case = TRUE))
  pos <- (n - 1L) %% GUIDES_PER_SERIES + 1L
  data.frame(
    Reference_ID = rid,
    series_idx = (n - 1L) %/% GUIDES_PER_SERIES + 1L,
    pos_in_series = pos,
    mutations_nt = pos * MUTATION_STEP_NT,
    stringsAsFactors = FALSE
  )
}

# Map each top50 series to its matched positive-control Sequence_name.
build_top50_pos_map <- function(seq_path) {
  va <- as.data.frame(readxl::read_excel(seq_path), stringsAsFactors = FALSE)
  pos <- va[va$Category == "PosCtrl", c("Sequence_name", "Chromosome", "HG38.start", "HG38.end")]
  top50 <- va[classify(va$Sequence_name) == "top50series",
              c("Sequence_name", "Chromosome", "HG38.start", "HG38.end")]
  meta <- parse_top50_ids(top50$Sequence_name)
  top50$series_idx <- meta$series_idx
  top50$pos_in_series <- meta$pos_in_series
  base <- top50[top50$pos_in_series == 1L, c("series_idx", "Chromosome", "HG38.start", "HG38.end")]
  m <- merge(base, pos, by = c("Chromosome", "HG38.start", "HG38.end"))
  out <- data.frame(
    series_idx = m$series_idx,
    pos_ctrl_id = m$Sequence_name,
    stringsAsFactors = FALSE
  )
  out[order(out$series_idx), , drop = FALSE]
}

series_with_pos_ctrl_anchor <- function(filt_path, top50_map) {
  filt <- read.table(filt_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  ctrl_ids <- filt$Reference_ID[classify(filt$Reference_ID) != "top50series"]
  top50_map$series_idx[top50_map$pos_ctrl_id %in% ctrl_ids]
}

mean_corrected_rate <- function(df, reference_id) {
  ref_rows <- df[as.character(df$Reference_ID) == as.character(reference_id), , drop = FALSE]
  if (nrow(ref_rows) == 0) return(NA_real_)
  avg_row <- ref_rows[as.character(ref_rows$replicate) == "avg", , drop = FALSE]
  if (nrow(avg_row) == 1L) return(avg_row$Corrected_TC_GA_Edit_Rate[1])
  rep_rows <- ref_rows[as.character(ref_rows$replicate) != "avg", , drop = FALSE]
  if (nrow(rep_rows) == 0) return(NA_real_)
  mean(rep_rows$Corrected_TC_GA_Edit_Rate, na.rm = TRUE)
}

# Adjacent paired Wilcoxon tests across mutation levels (BH-adjusted).
compute_pairwise_mutation_wilcox <- function(data, y_col) {
  mutation_points <- seq(0, MAX_MUTATIONS_NT - MUTATION_STEP_NT, by = MUTATION_STEP_NT)
  rows <- lapply(mutation_points, function(m1) {
    m2 <- m1 + MUTATION_STEP_NT
    d1 <- data[data$mutations_nt == m1, c("series_idx", y_col), drop = FALSE]
    names(d1)[2] <- "val_m1"
    d2 <- data[data$mutations_nt == m2, c("series_idx", y_col), drop = FALSE]
    names(d2)[2] <- "val_m2"
    paired <- merge(d1, d2, by = "series_idx")
    n_pairs <- nrow(paired)
    if (n_pairs < 3L) {
      return(data.frame(
        mutation_from_nt = m1, mutation_to_nt = m2,
        comparison = paste0(m1, " vs ", m2), n_pairs = n_pairs,
        wilcoxon_statistic = NA_real_, wilcoxon_p_value = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    wt <- suppressWarnings(stats::wilcox.test(paired$val_m2, paired$val_m1, paired = TRUE, exact = FALSE))
    data.frame(
      mutation_from_nt = m1, mutation_to_nt = m2,
      comparison = paste0(m1, " vs ", m2), n_pairs = n_pairs,
      wilcoxon_statistic = unname(wt$statistic), wilcoxon_p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$wilcoxon_p_adj_bh <- stats::p.adjust(out$wilcoxon_p_value, method = "BH")
  out$stars <- ifelse(
    is.na(out$wilcoxon_p_adj_bh), "",
    ifelse(out$wilcoxon_p_adj_bh < 0.001, "***",
           ifelse(out$wilcoxon_p_adj_bh < 0.01, "**",
                  ifelse(out$wilcoxon_p_adj_bh < 0.05, "*", "")))
  )
  out
}

# Trend (plasma series + black mean) and boxplot for mutation curves.
mutation_trend_box_plots <- function(df, y_col, ylab, hline0 = FALSE) {
  plot_df <- df
  plot_df$y <- plot_df[[y_col]]
  trend <- stats::aggregate(y ~ mutations_nt, data = plot_df, FUN = mean, na.rm = TRUE)
  names(trend)[2] <- "y_mean"
  trend <- trend[order(trend$mutations_nt), , drop = FALSE]

  hline <- if (hline0) {
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40")
  } else {
    NULL
  }

  p_trend <- ggplot2::ggplot(plot_df, ggplot2::aes(x = mutations_nt, y = y)) +
    hline +
    ggplot2::geom_line(ggplot2::aes(group = series_idx, color = series_idx),
                       alpha = 0.5, linewidth = 0.45) +
    ggplot2::scale_color_viridis_c(option = "plasma", guide = "none") +
    ggplot2::geom_line(data = trend, ggplot2::aes(x = mutations_nt, y = y_mean),
                       inherit.aes = FALSE, color = "black", linewidth = 1) +
    ggplot2::labs(x = "Mutated bases (nt)", y = ylab) +
    ggplot2::theme_classic(base_size = 12)

  p_box <- ggplot2::ggplot(plot_df, ggplot2::aes(x = factor(mutations_nt), y = y)) +
    hline +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.35, size = 0.8) +
    ggplot2::labs(x = "Mutated bases (nt)", y = ylab) +
    ggplot2::theme_classic(base_size = 12)

  list(trend = p_trend, box = p_box)
}
