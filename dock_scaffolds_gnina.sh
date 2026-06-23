#!/bin/bash
#===============================================================================
# GNINA CNN Docking of 15 MCS Scaffolds into GID1 (2ZSH)
#
# Docks each bare scaffold into the GA3 binding pocket using GNINA's CNN
# scoring to find optimal 3D poses for DiffGUI inpainting.
#
# Usage (on Expanse or any machine with GNINA + GPU):
#   bash dock_scaffolds_gnina.sh
#
# Or as a Slurm job:
#   sbatch dock_scaffolds_gnina.sh
#
#SBATCH --job-name=gnina_scaffolds
#SBATCH --account=iit135
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=scaffold_docking/gnina_dock_%j.out
#===============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DOCK_DIR="$SCRIPT_DIR/scaffold_docking"
RECEPTOR="$SCRIPT_DIR/2ZSH.pdb"
POCKET="$SCRIPT_DIR/2ZSH_pocket10.pdb"

# GA3 binding site center and box size
CENTER_X=51.033
CENTER_Y=59.452
CENTER_Z=37.370
SIZE_X=20
SIZE_Y=20
SIZE_Z=20

# Number of poses to keep per scaffold
NUM_POSES=3
EXHAUSTIVENESS=16

# Output directory for docked poses
RESULTS_DIR="$DOCK_DIR/docked_poses"
mkdir -p "$RESULTS_DIR"

echo "=== GNINA CNN Docking: 15 Scaffolds into GID1 (2ZSH) ==="
echo "Receptor: $RECEPTOR"
echo "Center: ($CENTER_X, $CENTER_Y, $CENTER_Z)"
echo "Box: ${SIZE_X}x${SIZE_Y}x${SIZE_Z} A"
echo "Poses per scaffold: $NUM_POSES"
echo ""

# Check for GNINA
if ! command -v gnina &> /dev/null; then
    echo "GNINA not found. Install or load module first."
    echo "On Expanse: module load gnina (or download from https://github.com/gnina/gnina/releases)"
    exit 1
fi

for SDF in "$DOCK_DIR"/*_3d.sdf; do
    BASENAME=$(basename "$SDF" _3d.sdf)
    OUT_SDF="$RESULTS_DIR/${BASENAME}_docked.sdf"
    OUT_LOG="$RESULTS_DIR/${BASENAME}_docking.log"

    echo "Docking: $BASENAME"

    gnina \
        --receptor "$RECEPTOR" \
        --ligand "$SDF" \
        --center_x $CENTER_X \
        --center_y $CENTER_Y \
        --center_z $CENTER_Z \
        --size_x $SIZE_X \
        --size_y $SIZE_Y \
        --size_z $SIZE_Z \
        --num_modes $NUM_POSES \
        --exhaustiveness $EXHAUSTIVENESS \
        --cnn_scoring rescore \
        --out "$OUT_SDF" \
        --log "$OUT_LOG" \
        2>&1 | tail -5

    echo "  -> $OUT_SDF"
    echo ""
done

echo "=== Docking Complete ==="
echo "Results in: $RESULTS_DIR/"
echo ""
echo "Next steps:"
echo "  1. Download the docked_poses/*.sdf files"
echo "  2. Open in PyMOL with 2ZSH.pdb to inspect poses"
echo "  3. Select best pose per scaffold for DiffGUI inpainting"

# Generate a combined SDF for easy PyMOL viewing
COMBINED="$RESULTS_DIR/all_scaffolds_docked.sdf"
cat "$RESULTS_DIR"/*_docked.sdf > "$COMBINED" 2>/dev/null
echo "Combined file: $COMBINED"
