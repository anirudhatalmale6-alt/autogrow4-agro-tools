#!/bin/bash
#=====================================================================
# Microenvironment Analysis - Task 2
# SDSC Expanse - Slurm Submission Script
#
# Sidechain centroid integration: void volume, electrostatic potential,
# SASA, and hydrophobicity for 14 key GID1 binding pocket residues.
#
# Usage:
#   bash submit_microenvironment.sh
#=====================================================================

PROJECT_DIR=~/final_project
OUTPUT_DIR="${PROJECT_DIR}/pocket_analysis"
SCRIPT="${PROJECT_DIR}/microenvironment_analysis.py"
RESULTS_DIR="${OUTPUT_DIR}/logs"

PARTITION="shared"
ACCOUNT="--account=iit135"
NTASKS=1
MEM="8G"
TIME="01:00:00"

mkdir -p "${OUTPUT_DIR}" "${RESULTS_DIR}"

echo "============================================"
echo "Microenvironment Analysis - Task 2"
echo "============================================"
echo "Output: ${OUTPUT_DIR}"
echo ""

# Pre-download PDB files on login node (compute nodes have no internet)
echo "Pre-downloading PDB files..."
for PDB in 2ZSH 2ZSI 3ED1 3EBL; do
    PDB_FILE="${OUTPUT_DIR}/${PDB}.pdb"
    if [ -f "${PDB_FILE}" ] && [ $(stat -c%s "${PDB_FILE}" 2>/dev/null || echo 0) -gt 1000 ]; then
        echo "  ${PDB}.pdb already exists, skipping"
    else
        echo "  Downloading ${PDB}..."
        wget -q "https://files.rcsb.org/download/${PDB}.pdb" -O "${PDB_FILE}"
        if [ $? -eq 0 ] && [ -s "${PDB_FILE}" ]; then
            echo "  ${PDB}.pdb downloaded OK"
        else
            echo "  ERROR: Failed to download ${PDB}.pdb"
            exit 1
        fi
    fi
done
echo "All PDB files ready."
echo ""

SLURM_SCRIPT="${OUTPUT_DIR}/slurm_microenvironment.sh"

cat > "${SLURM_SCRIPT}" << 'SLURM_EOF'
#!/bin/bash
#SBATCH --job-name=microenv
#SBATCH --partition=PARTITION_VAL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=NTASKS_VAL
#SBATCH --mem=MEM_VAL
#SBATCH --time=TIME_VAL
#SBATCH --output=OUTPUT_DIR_VAL/logs/microenv_%j.out
#SBATCH --error=OUTPUT_DIR_VAL/logs/microenv_%j.err
ACCOUNT_VAL

module purge
module load cpu/0.15.4
module load gcc/10.2.0

echo "============================================"
echo "Microenvironment Analysis - Task 2"
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
echo "Check output: tail -f ${RESULTS_DIR}/microenv_*.out"
echo "Results will be in: ${OUTPUT_DIR}/microenvironment_report.txt"
echo "Plots: microenvironment_deltav_deltae.png, microenvironment_landscape.png,"
echo "       microenvironment_summary.png, microenvironment_absolute.png"
