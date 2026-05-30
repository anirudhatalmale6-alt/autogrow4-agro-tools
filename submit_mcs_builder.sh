#!/bin/bash
#=====================================================================
# MCS Database Builder - Phase 2
# SDSC Expanse - Slurm Submission Script (CPU)
#
# Builds 15 scaffold-organized seed libraries from supplier catalog.
# CPU-only (no GPU needed - pure RDKit cheminformatics).
#
# Usage:
#   bash submit_mcs_builder.sh <input_file.smi>
#   bash submit_mcs_builder.sh <input_file.smi> quinoline  # single scaffold
#   bash submit_mcs_builder.sh --demo                      # test run
#=====================================================================

PROJECT_DIR=~/final_project
SCRIPT="${PROJECT_DIR}/mcs_database_builder.py"
OUTPUT_DIR="${PROJECT_DIR}/mcs_databases"
RESULTS_DIR="${OUTPUT_DIR}/logs"
ACCOUNT="--account=iit135"

INPUT_FILE="${1}"
SCAFFOLD="${2:-}"

if [ "$INPUT_FILE" = "--demo" ]; then
    MODE="demo"
    PYTHON_ARGS="--demo --output_dir ${OUTPUT_DIR}"
    MEM="8G"
    TIME="00:30:00"
elif [ -z "$INPUT_FILE" ]; then
    echo "Usage:"
    echo "  bash submit_mcs_builder.sh <input_file.smi>"
    echo "  bash submit_mcs_builder.sh <input_file.smi> quinoline"
    echo "  bash submit_mcs_builder.sh --demo"
    exit 1
else
    if [ ! -f "$INPUT_FILE" ]; then
        echo "ERROR: Input file not found: $INPUT_FILE"
        exit 1
    fi
    MODE="full"
    PYTHON_ARGS="--input ${INPUT_FILE} --output_dir ${OUTPUT_DIR}"
    if [ -n "$SCAFFOLD" ]; then
        PYTHON_ARGS="${PYTHON_ARGS} --scaffold ${SCAFFOLD}"
    fi
    MEM="32G"
    TIME="04:00:00"
fi

mkdir -p "${OUTPUT_DIR}" "${RESULTS_DIR}"

echo "============================================"
echo "MCS Database Builder - Phase 2"
echo "============================================"
echo "Mode: ${MODE}"
echo "Output: ${OUTPUT_DIR}"
if [ "$MODE" = "full" ]; then
    echo "Input: ${INPUT_FILE}"
    [ -n "$SCAFFOLD" ] && echo "Scaffold: ${SCAFFOLD}"
fi
echo ""

SLURM_SCRIPT="${OUTPUT_DIR}/slurm_mcs_builder.sh"

cat > "${SLURM_SCRIPT}" << 'SLURM_EOF'
#!/bin/bash
#SBATCH --job-name=mcsdb
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=MEM_VAL
#SBATCH --time=TIME_VAL
#SBATCH --output=OUTPUT_DIR_VAL/logs/mcsdb_%j.out
#SBATCH --error=OUTPUT_DIR_VAL/logs/mcsdb_%j.err
ACCOUNT_VAL

module purge
module load cpu/0.15.4
module load gcc/10.2.0

echo "============================================"
echo "MCS Database Builder - Phase 2"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "============================================"

source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/.conda/etc/profile.d/conda.sh 2>/dev/null
conda activate ag4_full 2>/dev/null

python3 -c "from rdkit import Chem" 2>/dev/null || { echo "ERROR: RDKit not found"; exit 1; }

python3 SCRIPT_VAL PYTHON_ARGS_VAL

echo ""
echo "============================================"
echo "Completed: $(date)"
echo "Results in: OUTPUT_DIR_VAL/"
echo "============================================"
SLURM_EOF

sed -i "s|MEM_VAL|${MEM}|g" "${SLURM_SCRIPT}"
sed -i "s|TIME_VAL|${TIME}|g" "${SLURM_SCRIPT}"
sed -i "s|ACCOUNT_VAL|#SBATCH ${ACCOUNT}|g" "${SLURM_SCRIPT}"
sed -i "s|OUTPUT_DIR_VAL|${OUTPUT_DIR}|g" "${SLURM_SCRIPT}"
sed -i "s|SCRIPT_VAL|${SCRIPT}|g" "${SLURM_SCRIPT}"
sed -i "s|PYTHON_ARGS_VAL|${PYTHON_ARGS}|g" "${SLURM_SCRIPT}"

echo "Submitting: sbatch ${SLURM_SCRIPT}"
sbatch "${SLURM_SCRIPT}"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Check output: tail -f ${RESULTS_DIR}/mcsdb_*.out"
echo "Results: ${OUTPUT_DIR}/pipeline_report.txt"
echo "         ${OUTPUT_DIR}/pipeline_results.json"
