#!/bin/bash
#=====================================================================
# Native Redocking Optimization - SDSC Expanse
# Finds optimal box dimensions for AutoDock Vina by minimizing RMSD
# between docked poses and crystallographic coordinates.
#
# Usage:
#   bash submit_redocking.sh          # Run both 2ZSH and 2ZSI
#   bash submit_redocking.sh 2ZSH     # Run only 2ZSH (GA3)
#   bash submit_redocking.sh 2ZSI     # Run only 2ZSI (GA4)
#=====================================================================

PROJECT_DIR=~/final_project
OUTPUT_DIR="${PROJECT_DIR}/redocking_optimization"
SCRIPT="${PROJECT_DIR}/redocking_optimizer.py"
RESULTS_DIR="${OUTPUT_DIR}/logs"

PDB_CHOICE="${1:-both}"
PARTITION="shared"
ACCOUNT="--account=iit135"
NTASKS=4
MEM="16G"
TIME="12:00:00"

mkdir -p "${RESULTS_DIR}"

SLURM_SCRIPT="${OUTPUT_DIR}/slurm_redocking.sh"

cat > "${SLURM_SCRIPT}" << 'SLURM_EOF'
#!/bin/bash
#SBATCH --job-name=redock_opt
#SBATCH --partition=PARTITION_VAL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=NTASKS_VAL
#SBATCH --mem=MEM_VAL
#SBATCH --time=TIME_VAL
#SBATCH --output=OUTPUT_DIR_VAL/logs/redocking_%j.out
#SBATCH --error=OUTPUT_DIR_VAL/logs/redocking_%j.err
ACCOUNT_VAL

module purge
module load cpu/0.15.4
module load gcc/10.2.0
module load openbabel/3.0.0

export PATH=~/final_project/autogrow4/autogrow/docking/docking_executables/vina/autodock_vina_1_1_2_linux_x86/bin:$PATH

echo "============================================"
echo "Native Redocking Box Optimization"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "============================================"

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/.conda/etc/profile.d/conda.sh 2>/dev/null
conda activate ag4_full 2>/dev/null

# Install pdb2pqr if not available
python3 -c "import pdb2pqr" 2>/dev/null || pip install pdb2pqr --quiet

python3 SCRIPT_VAL \
    --output-dir OUTPUT_DIR_VAL \
    --pdb PDB_VAL \
    --exhaustiveness 20

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
sed -i "s|PDB_VAL|${PDB_CHOICE}|g" "${SLURM_SCRIPT}"

echo "============================================"
echo "Native Redocking Box Optimization"
echo "============================================"
echo "PDB: ${PDB_CHOICE}"
echo "Output: ${OUTPUT_DIR}"
echo "Partition: ${PARTITION}"
echo "Resources: ${NTASKS} tasks, ${MEM} RAM, ${TIME}"
echo "============================================"
echo ""
echo "Submitting: sbatch ${SLURM_SCRIPT}"
sbatch "${SLURM_SCRIPT}"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Check output: tail -f ${RESULTS_DIR}/redocking_*.out"
echo "Results will be in: ${OUTPUT_DIR}/redocking_report.txt"
