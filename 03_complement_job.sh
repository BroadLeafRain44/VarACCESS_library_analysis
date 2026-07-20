#!/bin/bash
#SBATCH --job-name=ddda_complement
#SBATCH --account=def-glettre
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=03_complement_%j.out
#SBATCH --error=03_complement_%j.err

cd ~/links/projects/rrg-glettre/jordy2/programs/
source dna_env/bin/activate

module load StdEnv/2023
module load python/3.10

export DISABLE_ZLIB_NG=1

WORK_DIR=~/links/projects/rrg-glettre/jordy2/library/2026_02_09_library
OUTPUT_DIR=${WORK_DIR}/output
TMP_DIR=${SLURM_TMPDIR:-/tmp}

JOB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${JOB_DIR}/03_fastq_to_complement.py"
if [ ! -f "${PY_SCRIPT}" ]; then
    PY_SCRIPT="${WORK_DIR}/script/03_fastq_to_complement.py"
fi
if [ ! -f "${PY_SCRIPT}" ]; then
    echo "Error: 03_fastq_to_complement.py not found" >&2
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

MERGED_FASTQ=${OUTPUT_DIR}/NS.X0276.003.IndexR1_01.sample1.merged.fastq
TMP_COMPLEMENT="${TMP_DIR}/complement_${SLURM_JOB_ID}.fastq"
MERGED_COMPLEMENT=${OUTPUT_DIR}/NS.X0276.003.IndexR1_01.sample1.merged.complement.fastq

if [ ! -f "${MERGED_FASTQ}" ]; then
    echo "Error: Merged FASTQ file not found!"
    echo "Expected: ${MERGED_FASTQ}"
    exit 1
fi

echo "Create complement of merged FASTQ"
python3 "${PY_SCRIPT}" "${MERGED_FASTQ}" "${TMP_COMPLEMENT}"
mv "${TMP_COMPLEMENT}" "${MERGED_COMPLEMENT}"

echo "Complement creation complete"
echo "Output file: ${MERGED_COMPLEMENT}"
echo "Job completed at: $(date)"
