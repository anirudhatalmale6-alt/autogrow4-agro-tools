#!/bin/bash
#SBATCH --job-name=diffgui_frag
#SBATCH --account=iit135
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-15
#SBATCH --output=diffgui_run/logs/diffgui_%A_%a.out
#SBATCH --error=diffgui_run/logs/diffgui_%A_%a.err

module purge
module load cpu/0.17.3b
module load gpu/0.17.3b
module load anaconda3/2021.05
conda activate diffgui

DIFFGUI_DIR="/expanse/lustre/scratch/$USER/temp_project/DiffGui"
cd "$DIFFGUI_DIR"

SCAFFOLD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" diffgui_run/scaffold_names.txt)
CONFIG="diffgui_run/configs/${SCAFFOLD}_sample.yml"
OUTDIR="diffgui_run/outputs/${SCAFFOLD}"

echo "=== DiffGUI frag_cond: $SCAFFOLD ==="
echo "Config: $CONFIG"
echo "Output: $OUTDIR"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null)"
echo ""

mkdir -p "$OUTDIR"

python scripts/sample.py \
    --outdir "$OUTDIR" \
    --config "$CONFIG" \
    --device cuda:0

echo ""
echo "=== Done: $SCAFFOLD ==="
ls -la "$OUTDIR/" 2>/dev/null
