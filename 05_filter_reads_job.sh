#!/bin/bash
#SBATCH --job-name=ddda_filter_reads
#SBATCH --account=def-glettre
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=05_filter_reads_%j.out
#SBATCH --error=05_filter_reads_%j.err

cd ~/links/projects/rrg-glettre/jordy2/programs/
source dna_env/bin/activate

module load StdEnv/2023
module load python/3.10
module load samtools

export DISABLE_ZLIB_NG=1

WORK_DIR=~/links/projects/rrg-glettre/jordy2/library/2026_02_09_library
OUTPUT_DIR=${WORK_DIR}/output
TMP_DIR=${SLURM_TMPDIR:-/tmp}

JOB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${JOB_DIR}/05_filter_bam_by_read_length.py"
if [ ! -f "${PY_SCRIPT}" ]; then
    PY_SCRIPT="${WORK_DIR}/script/05_filter_bam_by_read_length.py"
fi
if [ ! -f "${PY_SCRIPT}" ]; then
    echo "Error: 05_filter_bam_by_read_length.py not found" >&2
    echo "Looked in: ${JOB_DIR} and ${WORK_DIR}/script" >&2
    exit 1
fi

mkdir -p ${OUTPUT_DIR}
cd ${WORK_DIR}

echo "Job started at: $(date)"
echo "Running on node: $SLURM_NODELIST"
echo "Job ID: $SLURM_JOB_ID"
echo "Using temporary directory: ${TMP_DIR}"
echo "Python script: ${PY_SCRIPT}"

BAM_INPUT=${OUTPUT_DIR}/NS.X0276.003.IndexR1_01.sample1.aligned.sorted.bam
BAM_FILTERED_TMP="${TMP_DIR}/filtered_${SLURM_JOB_ID}.bam"
BAM_FILTERED_SORTED_TMP="${TMP_DIR}/filtered_sorted_${SLURM_JOB_ID}.bam"
BAM_FILTERED=${OUTPUT_DIR}/NS.X0276.003.IndexR1_01.sample1.filtered_189_190bp.bam
BAM_FILTERED_SORTED=${OUTPUT_DIR}/NS.X0276.003.IndexR1_01.sample1.filtered_189_190bp.sorted.bam

TARGET_LENGTH_1=189
TARGET_LENGTH_2=190
COUNTS_FILE=${OUTPUT_DIR}/read_length_counts.txt

if [ ! -f "${BAM_INPUT}" ]; then
    echo "Error: Input BAM file not found!"
    echo "Expected: ${BAM_INPUT}"
    exit 1
fi

echo "Counting reads by length"
samtools view ${BAM_INPUT} | python3 "${PY_SCRIPT}" count > "${COUNTS_FILE}"
echo "Read length counts saved to: ${COUNTS_FILE}"
cat "${COUNTS_FILE}"

echo ""
echo "Filtering reads to keep only ${TARGET_LENGTH_1} bp and ${TARGET_LENGTH_2} bp reads"
samtools view -h ${BAM_INPUT} | \
    python3 "${PY_SCRIPT}" filter --lengths ${TARGET_LENGTH_1} ${TARGET_LENGTH_2} | \
    samtools view -bS - > "${BAM_FILTERED_TMP}"

mv "${BAM_FILTERED_TMP}" "${BAM_FILTERED}"

echo "Sorting filtered BAM file"
samtools sort "${BAM_FILTERED}" -o "${BAM_FILTERED_SORTED_TMP}"
samtools index "${BAM_FILTERED_SORTED_TMP}"

mv "${BAM_FILTERED_SORTED_TMP}" "${BAM_FILTERED_SORTED}"
if [ -f "${BAM_FILTERED_SORTED_TMP}.bai" ]; then
    mv "${BAM_FILTERED_SORTED_TMP}.bai" "${BAM_FILTERED_SORTED}.bai"
fi

echo "Filtered BAM file: ${BAM_FILTERED_SORTED}"
samtools flagstat ${BAM_FILTERED_SORTED}

echo "Read filtering complete"
echo "Output files:"
echo "${BAM_FILTERED}"
echo "${BAM_FILTERED_SORTED}"
echo "Job completed at: $(date)"
