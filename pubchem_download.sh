#!/bin/bash
#SBATCH --account=iit135
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --job-name="pubchem_dl"
#SBATCH --output="pubchem_dl.%j.%N.out"
#SBATCH --export=ALL

module purge
module load cpu

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ag4_full

cd ~/final_project

echo "Starting PubChem Agrochemical Library Builder"
echo "Date: $(date)"
echo ""

python pubchem_agro_library.py \
    --output autogrow4/source_compounds/agrochemical_library.smi \
    --skip-ghs \
    --resume

echo ""
echo "Finished: $(date)"
echo "Output file:"
ls -lh autogrow4/source_compounds/agrochemical_library.smi 2>/dev/null
wc -l autogrow4/source_compounds/agrochemical_library.smi 2>/dev/null
