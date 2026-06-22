#!/bin/bash
#SBATCH --job-name=diffgui_sample
#SBATCH --account=iit135
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=diffgui_%j.out
#SBATCH --error=diffgui_%j.err

#===============================================================================
# DiffGUI Sampling Job for SDSC Expanse
# Generates de novo ligands for a target protein pocket
#
# Usage:
#   sbatch submit_diffgui_sample.sh
#   sbatch submit_diffgui_sample.sh --export=POCKET=path/to/pocket.pdb,NUM_MOLS=200
#===============================================================================

module purge
module load cpu/0.17.3b
module load gpu/0.17.3b
module load anaconda3/2021.05

conda activate diffgui

WORK_DIR="/expanse/lustre/scratch/$USER/temp_project/DiffGui"
cd "$WORK_DIR"

# Allow override via environment variables
POCKET=${POCKET:-"sample/3ztx_pocket.pdb"}
NUM_MOLS=${NUM_MOLS:-10}
CONFIG=${CONFIG:-"configs/sample/sample.yml"}
OUTDIR=${OUTDIR:-"outputs"}
GUI_STRENGTH=${GUI_STRENGTH:-3.0}

echo "=== DiffGUI Sampling ==="
echo "Pocket: $POCKET"
echo "Num mols: $NUM_MOLS"
echo "Config: $CONFIG"
echo "Output: $OUTDIR"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'unknown')"
echo ""

python scripts/sample.py \
    --outdir "$OUTDIR" \
    --config "$CONFIG" \
    --device cuda:0

echo ""
echo "=== Sampling Complete ==="
echo "Results in: $WORK_DIR/$OUTDIR/"
ls -la "$WORK_DIR/$OUTDIR/" 2>/dev/null || echo "No output files found"
