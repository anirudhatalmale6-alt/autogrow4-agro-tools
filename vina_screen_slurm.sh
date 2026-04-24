#!/bin/bash
#SBATCH --job-name=vina_screen
#SBATCH --account=iit135
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=vina_screen.%j.out
#SBATCH --error=vina_screen.%j.err
#SBATCH --requeue

module purge
module load cpu

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ag4_full

cd ~/final_project

# ── Configuration ──
# Change INPUT_SMI to screen different databases
INPUT_SMI="autogrow4/source_compounds/agrochemical_library_v2_binD_200to250.smi"
RECEPTOR="receptors/2ZSH.pdbqt"
OUTPUT_PREFIX="screen_binD"

# Binding box (centered on GA3 ligand in GID1)
CENTER_X=51.0
CENTER_Y=59.5
CENTER_Z=37.4
SIZE_X=20.0
SIZE_Y=20.0
SIZE_Z=20.0

# Screening parameters
EXHAUSTIVENESS=1        # Low for fast screening (increase for accuracy)
TOP_PERCENT=10          # Keep top 10%
NPROCS=${SLURM_CPUS_PER_TASK:-100}
TIMEOUT=120             # Per-compound timeout (seconds)

echo "Virtual Screen Configuration"
echo "  Input: ${INPUT_SMI}"
echo "  Receptor: ${RECEPTOR}"
echo "  Workers: ${NPROCS}"
echo "  Date: $(date)"
echo ""

python -u vina_screen.py \
    --input "${INPUT_SMI}" \
    --receptor "${RECEPTOR}" \
    --output "${OUTPUT_PREFIX}" \
    --center_x ${CENTER_X} \
    --center_y ${CENTER_Y} \
    --center_z ${CENTER_Z} \
    --size_x ${SIZE_X} \
    --size_y ${SIZE_Y} \
    --size_z ${SIZE_Z} \
    --exhaustiveness ${EXHAUSTIVENESS} \
    --top_percent ${TOP_PERCENT} \
    --nprocs ${NPROCS} \
    --timeout ${TIMEOUT} \
    --resume

echo ""
echo "Finished: $(date)"
echo "Output files:"
ls -lh ${OUTPUT_PREFIX}* 2>/dev/null
