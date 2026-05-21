#!/bin/bash
#=====================================================================
# Cross-Docking IFP Matrix - Task 3
# SDSC Expanse - Slurm Submission Script (GPU)
#
# Rigid cross-docking validation: docks GA3 and GA4 across all 4 GID1
# structures using GNINA CNN scoring, computes interaction fingerprints,
# and calculates Tanimoto similarity matrix.
#
# Usage:
#   bash submit_cross_docking.sh          # GPU (recommended)
#   bash submit_cross_docking.sh cpu      # CPU fallback
#=====================================================================

PROJECT_DIR=~/final_project
OUTPUT_DIR="${PROJECT_DIR}/cross_docking"
POCKET_DIR="${PROJECT_DIR}/pocket_analysis"
SCRIPT="${PROJECT_DIR}/cross_docking_ifp.py"
RESULTS_DIR="${OUTPUT_DIR}/logs"

USE_GPU="${1:-gpu}"
ACCOUNT="--account=iit135"
MEM="16G"
TIME="04:00:00"

if [ "$USE_GPU" = "cpu" ]; then
    PARTITION="shared"
    NTASKS=4
    GPU_LINE=""
    TIME="08:00:00"
    echo "Mode: CPU (slower, no GPU acceleration)"
else
    PARTITION="gpu-shared"
    NTASKS=4
    GPU_LINE="#SBATCH --gpus=1"
    echo "Mode: GPU (recommended for GNINA CNN scoring)"
fi

mkdir -p "${OUTPUT_DIR}" "${RESULTS_DIR}" "${POCKET_DIR}"

echo "============================================"
echo "Cross-Docking IFP Matrix - Task 3"
echo "============================================"
echo "Partition: ${PARTITION}"
echo "Output: ${OUTPUT_DIR}"
echo ""

# Pre-download PDB files on login node
echo "Pre-downloading PDB files..."
for PDB in 2ZSH 2ZSI 3ED1 3EBL; do
    PDB_FILE="${POCKET_DIR}/${PDB}.pdb"
    if [ -f "${PDB_FILE}" ] && [ $(stat -c%s "${PDB_FILE}" 2>/dev/null || echo 0) -gt 1000 ]; then
        echo "  ${PDB}.pdb already exists"
    else
        echo "  Downloading ${PDB}..."
        wget -q "https://files.rcsb.org/download/${PDB}.pdb" -O "${PDB_FILE}"
        if [ $? -ne 0 ] || [ ! -s "${PDB_FILE}" ]; then
            echo "  ERROR: Failed to download ${PDB}.pdb"
            exit 1
        fi
    fi
done

# Pre-download GNINA if not present
GNINA_PATH="${PROJECT_DIR}/gnina"
if [ ! -x "$GNINA_PATH" ]; then
    echo "Downloading GNINA..."
    wget -q https://github.com/gnina/gnina/releases/download/v1.1/gnina -O "$GNINA_PATH"
    chmod +x "$GNINA_PATH"
    echo "GNINA downloaded."
else
    echo "GNINA already present."
fi
echo ""

SLURM_SCRIPT="${OUTPUT_DIR}/slurm_cross_docking.sh"

cat > "${SLURM_SCRIPT}" << 'SLURM_EOF'
#!/bin/bash
#SBATCH --job-name=crossdock
#SBATCH --partition=PARTITION_VAL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=NTASKS_VAL
#SBATCH --mem=MEM_VAL
#SBATCH --time=TIME_VAL
#SBATCH --output=OUTPUT_DIR_VAL/logs/crossdock_%j.out
#SBATCH --error=OUTPUT_DIR_VAL/logs/crossdock_%j.err
ACCOUNT_VAL
GPU_LINE_VAL

module purge
module load cpu/0.15.4
module load gcc/10.2.0

echo "============================================"
echo "Cross-Docking IFP Matrix - Task 3"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "============================================"

source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/.conda/etc/profile.d/conda.sh 2>/dev/null
conda activate ag4_full 2>/dev/null

python3 -c "import Bio" 2>/dev/null || pip install biopython --quiet
python3 -c "import matplotlib" 2>/dev/null || pip install matplotlib --quiet

echo "GNINA version:"
GNINA_PATH_VAL --version 2>&1 || echo "  (version check skipped)"
echo ""

python3 SCRIPT_VAL

echo ""
echo "============================================"
echo "Completed: $(date)"
echo "Results in: OUTPUT_DIR_VAL/"
echo "============================================"
SLURM_EOF

sed -i "s|PARTITION_VAL|${PARTITION}|g" "${SLURM_SCRIPT}"
sed -i "s|NTASKS_VAL|${NTASKS}|g" "${SLURM_SCRIPT}"
sed -i "s|MEM_VAL|${MEM}|g" "${SLURM_SCRIPT}"
sed -i "s|TIME_VAL|${TIME}|g" "${SLURM_SCRIPT}"
sed -i "s|ACCOUNT_VAL|#SBATCH ${ACCOUNT}|g" "${SLURM_SCRIPT}"
sed -i "s|OUTPUT_DIR_VAL|${OUTPUT_DIR}|g" "${SLURM_SCRIPT}"
sed -i "s|SCRIPT_VAL|${SCRIPT}|g" "${SLURM_SCRIPT}"
sed -i "s|GNINA_PATH_VAL|${GNINA_PATH}|g" "${SLURM_SCRIPT}"
sed -i "s|GPU_LINE_VAL|${GPU_LINE}|g" "${SLURM_SCRIPT}"

echo "Submitting: sbatch ${SLURM_SCRIPT}"
sbatch "${SLURM_SCRIPT}"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Check output: tail -f ${RESULTS_DIR}/crossdock_*.out"
echo "Results: ${OUTPUT_DIR}/cross_docking_ifp_report.txt"
echo "         ${OUTPUT_DIR}/cross_docking_results.json"
echo "Plots:   ${OUTPUT_DIR}/cross_docking_ifp_matrix.png"
