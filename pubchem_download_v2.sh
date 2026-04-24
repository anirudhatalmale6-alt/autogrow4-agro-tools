#!/bin/bash
#SBATCH --job-name=pubchem_agro_v2
#SBATCH --account=iit135
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=pubchem_dl_v2.%j.out
#SBATCH --error=pubchem_dl_v2.%j.err

module purge
module load cpu
source activate ag4_full

cd ~/final_project

python pubchem_agro_library_v2.py \
    --output autogrow4/source_compounds/agrochemical_library_v2.smi \
    --skip-ghs \
    --resume
