#!/bin/bash
#SBATCH --job-name=ddda_motifs_simple
#SBATCH --account=def-glettre
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=08_compare_reads_to_reference_with_motifs_%j.out
#SBATCH --error=08_compare_reads_to_reference_with_motifs_%j.err

module load StdEnv/2023
module load python/3.10

cd ~/links/projects/rrg-glettre/jordy2/programs/
source dna_env/bin/activate

export DISABLE_ZLIB_NG=1

cd ~/links/projects/rrg-glettre/jordy2/library/2026_02_09_library/

echo "Job started at: $(date)"
python3 -c "exec(open('script/08_compare_reads_to_reference_with_motifs.py').read())"
echo "Job completed at: $(date)"
