# VarACCESS library analysis

Scripts for processing and analyzing the **VarACCESS** oligo library: a SLURM sequencing pipeline (trim → merge → align → filter → motif comparison) plus R scripts for QC, edit-rate correction, and figure generation.

> **Note:** Personal machine/HPC paths have been replaced with placeholders (`/PATH/TO/YOUR/PROJECT`, `/PATH/TO/YOUR/PROGRAMS/`, `YOUR_ACCOUNT`, etc.). Search for `EDIT ME` and fill in your own paths, account, sample filenames, and references before running.

---

## Contents

### Sequencing pipeline (HPC / SLURM)

Run roughly in numerical order:

| Step | Files | What it does |
|------|--------|----------------|
| 1 | `01_cutadapt_job.sh` | Trim adapters from paired FASTQs |
| 2 | `02_flash_merge_job.sh` | Merge overlapping R1/R2 reads |
| 3 | `03_complement_job.sh`, `03_fastq_to_complement.py` | Reverse-complement reads as needed |
| 4 | `04_align_bam_job.sh` | Align to reference (BWA) |
| 5 | `05_filter_reads_job.sh`, `05_filter_bam_by_read_length.py` | Filter BAM by read length |
| 8 | `08_compare_reads_to_reference_with_motifs_job.sh`, `08_compare_reads_to_reference_with_motifs.py` | Compare reads to references / motifs |

### R analysis & plots

Shared helpers live in `library_helpers.R` (IDs, paths, top50series parsing, Excel/TSV readers).

| Script | Purpose |
|--------|---------|
| `filter_raw_data_read50_min2replicates.R` | Filter raw counts (min read depth / replicates) |
| `log_correct_GAM_TC_GA_edit_rate.R` | GAM log-correction of TC/GA edit rates |
| `log_cor_pos_neg_roc_simple.R` | ROC/AUC for positive vs negative controls |
| `plot_*.R` | Figures (edit rates, AlphaG, ChromBPNet, controls, top50series) |
| `series_resolvability_pooled_noise_raw_lfc.R` | Series resolvability / noise analysis |

---

## Requirements

**Pipeline (typical):** Python 3, cutadapt, FLASH, BWA, samtools, pysam; SLURM job scripts written for a Compute Canada–style environment.

**R packages used across scripts:** `ggplot2`, `dplyr`, `readxl`, `readr`, `patchwork`, `mgcv` (plus base R).

---

## How to use

1. Clone or download this repo.
2. Set placeholders marked `# EDIT ME`:
   - **R scripts:** `base_dir <- "/PATH/TO/YOUR/PROJECT"` (root with `data/`, `clean_outputs/`, etc.) and the `source(.../library_helpers.R)` path if needed.
   - **SLURM jobs:** `YOUR_ACCOUNT`, `/PATH/TO/YOUR/PROJECT`, `/PATH/TO/YOUR/PROGRAMS/`, sample FASTQ/BAM names, and reference FASTA names.
3. For pipeline jobs: submit with `sbatch` on your cluster (after adjusting account, modules, and file paths).
4. For R scripts: from a terminal, or in RStudio:

```bash
Rscript name_of_script.R
```

Most R scripts source `library_helpers.R` from `base_dir/git_scripts/` — keep that layout, or edit the `source(...)` line at the top of each script.

---

## License / reuse

Feel free to adapt these scripts for related library analyses. If you reuse them in a publication, a citation or acknowledgment is appreciated.
