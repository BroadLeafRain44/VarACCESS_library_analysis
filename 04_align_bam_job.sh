#!/bin/bash
#SBATCH --job-name=ddda_align
#SBATCH --account=def-glettre
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=04_align_%j.out
#SBATCH --error=04_align_%j.err

# Setup environment
cd ~/links/projects/rrg-glettre/jordy2/programs/
source dna_env/bin/activate

module load StdEnv/2023
module load python/3.10 
module load bwa 
module load samtools

# Disable zlib_ng to avoid "Bad address" errors on HPC filesystems
export DISABLE_ZLIB_NG=1

# Set paths
WORK_DIR=~/links/projects/rrg-glettre/jordy2/library/2026_02_09_library
REF_DIR=${WORK_DIR}/reference
OUTPUT_DIR=${WORK_DIR}/output

mkdir -p ${OUTPUT_DIR}

cd ${WORK_DIR}

# Use local scratch to avoid Lustre filesystem write issues
TMP_DIR=${SLURM_TMPDIR:-/tmp}

echo "Job started at: $(date)"
echo "Running on node: $SLURM_NODELIST"
echo "Job ID: $SLURM_JOB_ID"
echo "Using temporary directory: ${TMP_DIR}"

# Input FASTQ files (merged and complement from previous steps)
MERGED_FASTQ=${OUTPUT_DIR}/NS.X0276.003.IndexR1_01.sample1.merged.fastq
MERGED_COMPLEMENT=${OUTPUT_DIR}/NS.X0276.003.IndexR1_01.sample1.merged.complement.fastq

# Check if input files exist
if [ ! -f "${MERGED_FASTQ}" ]; then
    echo "Error: Merged FASTQ file not found!"
    echo "Expected: ${MERGED_FASTQ}"
    exit 1
fi

if [ ! -f "${MERGED_COMPLEMENT}" ]; then
    echo "Error: Complement FASTQ file not found!"
    echo "Expected: ${MERGED_COMPLEMENT}"
    exit 1
fi

# Reference files
REF_FASTA=${REF_DIR}/2026_02_09_reference.fasta
REF_COMPLEMENT=${REF_DIR}/2026_02_09_reference_complement.fasta

# Temporary file paths for main alignment
TMP_SAM="${TMP_DIR}/aligned_${SLURM_JOB_ID}.sam"
TMP_BAM="${TMP_DIR}/aligned_${SLURM_JOB_ID}.bam"
TMP_BAM_SORTED="${TMP_DIR}/aligned_sorted_${SLURM_JOB_ID}.bam"

# Temporary file paths for complement alignment
TMP_SAM_COMPLEMENT="${TMP_DIR}/aligned_complement_${SLURM_JOB_ID}.sam"
TMP_BAM_COMPLEMENT="${TMP_DIR}/aligned_complement_${SLURM_JOB_ID}.bam"
TMP_BAM_SORTED_COMPLEMENT="${TMP_DIR}/aligned_sorted_complement_${SLURM_JOB_ID}.bam"

# Final output file paths for main alignment
SAM_FILE=${OUTPUT_DIR}/NS.X0276.003.IndexR1_01.sample1.aligned.sam
BAM_FILE=${OUTPUT_DIR}/NS.X0276.003.IndexR1_01.sample1.aligned.bam
BAM_SORTED=${OUTPUT_DIR}/NS.X0276.003.IndexR1_01.sample1.aligned.sorted.bam

# Final output file paths for complement alignment
SAM_FILE_COMPLEMENT=${OUTPUT_DIR}/NS.X0276.003.IndexR1_01.sample1.aligned.complement.sam
BAM_FILE_COMPLEMENT=${OUTPUT_DIR}/NS.X0276.003.IndexR1_01.sample1.aligned.complement.bam
BAM_SORTED_COMPLEMENT=${OUTPUT_DIR}/NS.X0276.003.IndexR1_01.sample1.aligned.complement.sorted.bam

echo "Index reference sequences"
bwameth.py index ${REF_FASTA}
bwameth.py index ${REF_COMPLEMENT}

echo "Align merged reads to reference"
bwameth.py --reference ${REF_FASTA} ${MERGED_FASTQ} > "${TMP_SAM}"
mv "${TMP_SAM}" "${SAM_FILE}"

echo "Convert SAM to BAM and sort (main alignment)"
samtools view -bS ${SAM_FILE} > "${TMP_BAM}"
mv "${TMP_BAM}" "${BAM_FILE}"
samtools sort "${BAM_FILE}" -o "${TMP_BAM_SORTED}"
samtools index "${TMP_BAM_SORTED}"

# Move sorted BAM and index to final location
mv "${TMP_BAM_SORTED}" "${BAM_SORTED}"
if [ -f "${TMP_BAM_SORTED}.bai" ]; then
    mv "${TMP_BAM_SORTED}.bai" "${BAM_SORTED}.bai"
fi

echo "Align complement reads to complement reference"
bwameth.py --reference ${REF_COMPLEMENT} ${MERGED_COMPLEMENT} > "${TMP_SAM_COMPLEMENT}"
mv "${TMP_SAM_COMPLEMENT}" "${SAM_FILE_COMPLEMENT}"

echo "Convert SAM to BAM and sort (complement alignment)"
samtools view -bS ${SAM_FILE_COMPLEMENT} > "${TMP_BAM_COMPLEMENT}"
mv "${TMP_BAM_COMPLEMENT}" "${BAM_FILE_COMPLEMENT}"
samtools sort "${BAM_FILE_COMPLEMENT}" -o "${TMP_BAM_SORTED_COMPLEMENT}"
samtools index "${TMP_BAM_SORTED_COMPLEMENT}"

# Move sorted BAM and index to final location
mv "${TMP_BAM_SORTED_COMPLEMENT}" "${BAM_SORTED_COMPLEMENT}"
if [ -f "${TMP_BAM_SORTED_COMPLEMENT}.bai" ]; then
    mv "${TMP_BAM_SORTED_COMPLEMENT}.bai" "${BAM_SORTED_COMPLEMENT}.bai"
fi

echo "Check alignment statistics"
echo "Main alignment statistics:"
samtools flagstat ${BAM_SORTED}
echo ""
echo "Complement alignment statistics:"
samtools flagstat ${BAM_SORTED_COMPLEMENT}

echo "Alignment and BAM conversion complete"
echo "Main alignment output files:"
echo "${SAM_FILE}"
echo "${BAM_FILE}"
echo "${BAM_SORTED}"
echo ""
echo "Complement alignment output files:"
echo "${SAM_FILE_COMPLEMENT}"
echo "${BAM_FILE_COMPLEMENT}"
echo "${BAM_SORTED_COMPLEMENT}"
echo "Job completed at: $(date)"