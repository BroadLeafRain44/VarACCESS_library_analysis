#!/bin/bash
#SBATCH --job-name=ddda_flash
#SBATCH --account=YOUR_ACCOUNT  # EDIT ME
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=02_flash_%j.out
#SBATCH --error=02_flash_%j.err

# Setup environment
cd /PATH/TO/YOUR/PROGRAMS/  # EDIT ME
source dna_env/bin/activate  # EDIT ME: your venv

module load StdEnv/2023

# Disable zlib_ng to avoid "Bad address" errors on HPC filesystems
export DISABLE_ZLIB_NG=1

# Set paths
WORK_DIR=/PATH/TO/YOUR/PROJECT  # EDIT ME
FASTQ_DIR=${WORK_DIR}/fastq
OUTPUT_DIR=${WORK_DIR}/output

mkdir -p ${OUTPUT_DIR}

cd ${WORK_DIR}

# Use local scratch to avoid Lustre filesystem write issues
TMP_DIR=${SLURM_TMPDIR:-/tmp}

echo "Job started at: $(date)"
echo "Running on node: $SLURM_NODELIST"
echo "Job ID: $SLURM_JOB_ID"
echo "Using temporary directory: ${TMP_DIR}"

# Input FASTQ files (trimmed R1/R2 from previous step)
R1_TRIMMED=${FASTQ_DIR}/sample1_R1.trimmed.fastq
R2_TRIMMED=${FASTQ_DIR}/sample1_R2.trimmed.fastq

# Check if input files exist
if [ ! -f "${R1_TRIMMED}" ] || [ ! -f "${R2_TRIMMED}" ]; then
    echo "Error: Trimmed FASTQ files not found!"
    echo "Expected: ${R1_TRIMMED}"
    echo "Expected: ${R2_TRIMMED}"
    exit 1
fi

# Final output file path
MERGED_FASTQ=${OUTPUT_DIR}/sample1.merged.fastq

# FLASH binary
FLASH_BIN=/PATH/TO/YOUR/PROGRAMS/FLASH-1.2.11/flash

echo "Merge R1 and R2 reads using FLASH"
${FLASH_BIN} -t ${SLURM_CPUS_PER_TASK} -m 10 -M 150 \
  -d ${TMP_DIR} \
  -o merged_${SLURM_JOB_ID} \
  ${R1_TRIMMED} ${R2_TRIMMED}

# FLASH creates files with .extendedFrags.fastq extension
FLASH_TMP_OUTPUT=${TMP_DIR}/merged_${SLURM_JOB_ID}.extendedFrags.fastq

if [ ! -f "${FLASH_TMP_OUTPUT}" ]; then
    echo "Error: FLASH output file not found!"
    exit 1
fi

mv "${FLASH_TMP_OUTPUT}" "${MERGED_FASTQ}"

echo "FLASH merge complete"
echo "Output file: ${MERGED_FASTQ}"
echo "Job completed at: $(date)"
