#!/bin/bash
#=====================================================================
# Pocket RMSD Analysis - Universal Ligand Feasibility (Task 1)
# SDSC Expanse - Slurm Submission Script
#
# Aligns 4 GID1 crystal structures on 14 key binding pocket residues
# and calculates pairwise RMSD to determine if a universal ligand
# can bind both monocot (rice) and dicot (Arabidopsis) GID1.
#
# Usage:
#   bash submit_pocket_rmsd.sh
#=====================================================================

PROJECT_DIR=~/final_project
OUTPUT_DIR="${PROJECT_DIR}/pocket_analysis"
SCRIPT="${PROJECT_DIR}/pocket_rmsd.py"
RESULTS_DIR="${OUTPUT_DIR}/logs"

PARTITION="shared"
ACCOUNT="--account=iit135"
NTASKS=1
MEM="4G"
TIME="00:30:00"

mkdir -p "${OUTPUT_DIR}" "${RESULTS_DIR}"

echo "============================================"
echo "Pocket RMSD Analysis - Task 1"
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

SLURM_SCRIPT="${OUTPUT_DIR}/slurm_pocket_rmsd.sh"

cat > "${SLURM_SCRIPT}" << 'SLURM_EOF'
#!/bin/bash
#SBATCH --job-name=pocket_rmsd
#SBATCH --partition=PARTITION_VAL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=NTASKS_VAL
#SBATCH --mem=MEM_VAL
#SBATCH --time=TIME_VAL
#SBATCH --output=OUTPUT_DIR_VAL/logs/pocket_rmsd_%j.out
#SBATCH --error=OUTPUT_DIR_VAL/logs/pocket_rmsd_%j.err
ACCOUNT_VAL

module purge
module load cpu/0.15.4
module load gcc/10.2.0

echo "============================================"
echo "Pocket RMSD Analysis - Task 1"
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
echo "Check output: tail -f ${RESULTS_DIR}/pocket_rmsd_*.out"
echo "Results will be in: ${OUTPUT_DIR}/pocket_rmsd_report.txt"
echo "Plots will be in:   ${OUTPUT_DIR}/pocket_rmsd_analysis.png"
echo "                     ${OUTPUT_DIR}/pocket_rmsd_summary.png"
