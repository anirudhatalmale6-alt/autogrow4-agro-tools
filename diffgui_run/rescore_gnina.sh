#!/bin/bash
#SBATCH --job-name=gnina_rescore
#SBATCH --account=iit135
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=diffgui_run/logs/gnina_rescore_%j.out
#SBATCH --error=diffgui_run/logs/gnina_rescore_%j.err

# GNINA CNN + Vina rescoring of DiffGUI outputs against 2ZSH
# Requires gnina installed in the conda environment or available as module

module purge
module load cpu/0.17.3b
module load gpu/0.17.3b
module load anaconda3/2021.05
conda activate diffgui

OUTDIR="diffgui_run"
RECEPTOR="2ZSH.pdb"

echo "=== GNINA CNN Rescoring ==="
echo "Receptor: $RECEPTOR"
echo ""

RESCORE_DIR="$OUTDIR/gnina_rescores"
mkdir -p "$RESCORE_DIR"

for SCAFFOLD_DIR in "$OUTDIR"/outputs/*/; do
    SCAFFOLD=$(basename "$SCAFFOLD_DIR")
    echo "Rescoring: $SCAFFOLD"

    # Find all generated SDF files in the scaffold output
    for SDF in "$SCAFFOLD_DIR"/*.sdf; do
        [ -f "$SDF" ] || continue
        BASENAME=$(basename "$SDF" .sdf)
        OUT_SDF="$RESCORE_DIR/${SCAFFOLD}_${BASENAME}_rescored.sdf"

        gnina \
            -r "$RECEPTOR" \
            -l "$SDF" \
            --score_only \
            --cnn_scoring \
            -o "$OUT_SDF" \
            --seed 42 \
            2>/dev/null

        if [ -f "$OUT_SDF" ]; then
            echo "  $BASENAME: done"
        else
            echo "  $BASENAME: FAILED"
        fi
    done
done

echo ""
echo "=== Rescoring Complete ==="
echo "Results in: $RESCORE_DIR/"
ls "$RESCORE_DIR/"*.sdf 2>/dev/null | wc -l
echo "SDF files generated"
