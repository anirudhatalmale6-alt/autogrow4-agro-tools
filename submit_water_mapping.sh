#!/bin/bash
#=====================================================================
# Conserved Water Mapping - Task 5
# SDSC Expanse - Slurm Submission Script
#
# Maps crystallographic waters in the GID1 binding pocket, identifies
# conserved positions across 4 structures, and analyzes 4 specific
# water-mediated hydrogen bonding networks.
#
# Usage:
#   bash submit_water_mapping.sh
#=====================================================================

PROJECT_DIR=~/final_project
OUTPUT_DIR="${PROJECT_DIR}/water_mapping"
POCKET_DIR="${PROJECT_DIR}/pocket_analysis"
SCRIPT="${PROJECT_DIR}/conserved_water_mapping.py"
RESULTS_DIR="${OUTPUT_DIR}/logs"

PARTITION="shared"
ACCOUNT="--account=iit135"
NTASKS=1
MEM="8G"
TIME="01:00:00"

mkdir -p "${OUTPUT_DIR}" "${RESULTS_DIR}" "${POCKET_DIR}"

echo "============================================"
echo "Conserved Water Mapping - Task 5"
echo "============================================"
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
echo "All PDB files ready."
echo ""

SLURM_SCRIPT="${OUTPUT_DIR}/slurm_water_mapping.sh"

cat > "${SLURM_SCRIPT}" << 'SLURM_EOF'
#!/bin/bash
#SBATCH --job-name=watermap
#SBATCH --partition=PARTITION_VAL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=NTASKS_VAL
#SBATCH --mem=MEM_VAL
#SBATCH --time=TIME_VAL
#SBATCH --output=OUTPUT_DIR_VAL/logs/watermap_%j.out
#SBATCH --error=OUTPUT_DIR_VAL/logs/watermap_%j.err
ACCOUNT_VAL

module purge
module load cpu/0.15.4
module load gcc/10.2.0

echo "============================================"
echo "Conserved Water Mapping - Task 5"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "============================================"

source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/.conda/etc/profile.d/conda.sh 2>/dev/null
conda activate ag4_full 2>/dev/null

python3 -c "import Bio" 2>/dev/null || pip install biopython --quiet
python3 -c "import matplotlib" 2>/dev/null || pip install matplotlib --quiet

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

echo "Submitting: sbatch ${SLURM_SCRIPT}"
sbatch "${SLURM_SCRIPT}"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Check output: tail -f ${RESULTS_DIR}/watermap_*.out"
echo "Results: ${OUTPUT_DIR}/water_mapping_report.txt"
echo "         ${OUTPUT_DIR}/water_mapping_results.json"
echo "Plots:   ${OUTPUT_DIR}/water_mapping_analysis.png"
echo "         ${OUTPUT_DIR}/water_bridge_summary.png"
echo "Viz:     ${OUTPUT_DIR}/visualize_waters.py"
