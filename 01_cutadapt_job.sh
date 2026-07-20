#!/bin/bash
#SBATCH --job-name=ddda_cutadapt
#SBATCH --account=YOUR_ACCOUNT  # EDIT ME
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=01_cutadapt_%j.out
#SBATCH --error=01_cutadapt_%j.err

# Setup environment
cd /PATH/TO/YOUR/PROGRAMS/  # EDIT ME
source dna_env/bin/activate  # EDIT ME: your venv

module load StdEnv/2023
module load python/3.10

# Disable zlib_ng to avoid "Bad address" errors on HPC filesystems
export DISABLE_ZLIB_NG=1

# Set paths
WORK_DIR=/PATH/TO/YOUR/PROJECT  # EDIT ME
FASTQ_DIR=${WORK_DIR}/fastq

mkdir -p ${FASTQ_DIR}

cd ${WORK_DIR}

# Use local scratch to avoid Lustre filesystem write issues
TMP_DIR=${SLURM_TMPDIR:-/tmp}

echo "Job started at: $(date)"
echo "Running on node: $SLURM_NODELIST"
echo "Job ID: $SLURM_JOB_ID"
echo "Using temporary directory: ${TMP_DIR}"

# Input FASTQ files (raw R1/R2 from sequencer)
R1_INPUT=${FASTQ_DIR}/sample1_R1.fastq.gz
R2_INPUT=${FASTQ_DIR}/sample1_R2.fastq.gz

# Temporary file paths
TMP_R1_TRIMMED="${TMP_DIR}/R1_trimmed_${SLURM_JOB_ID}.fastq"
TMP_R2_TRIMMED="${TMP_DIR}/R2_trimmed_${SLURM_JOB_ID}.fastq"

# Final output file paths
R1_TRIMMED=${FASTQ_DIR}/sample1_R1.trimmed.fastq
R2_TRIMMED=${FASTQ_DIR}/sample1_R2.trimmed.fastq

echo "Trim adapters from reads"
cutadapt \
  -g "TTRTAYYTTAYYARYYAYYA" \
  -G "TRRTAAAAYRRYATARRRRT" \
  -m 30 \
  -q 20 \
  --quality-base=33 \
  --pair-filter=both \
  -o "${TMP_R1_TRIMMED}" \
  -p "${TMP_R2_TRIMMED}" \
  ${R1_INPUT} ${R2_INPUT}

# Move trimmed files to final location
echo "Moving trimmed files to final destination..."
mv "${TMP_R1_TRIMMED}" "${R1_TRIMMED}"
mv "${TMP_R2_TRIMMED}" "${R2_TRIMMED}"

echo "Cutadapt complete"
echo "Output files:"
echo "${R1_TRIMMED}"
echo "${R2_TRIMMED}"
echo "Job completed at: $(date)"
