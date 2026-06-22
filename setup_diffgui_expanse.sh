#!/bin/bash
#===============================================================================
# DiffGUI Setup Script for SDSC Expanse
# Target-aware 3D Molecular Generation (SE(3) Equivariant Diffusion Model)
# https://github.com/QiaoyuHu89/DiffGui
#
# Run this interactively on Expanse login node:
#   bash setup_diffgui_expanse.sh
#
# After setup, submit jobs with:
#   sbatch submit_diffgui_sample.sh
#===============================================================================

set -e

WORK_DIR="/expanse/lustre/scratch/$USER/temp_project"
DIFFGUI_DIR="$WORK_DIR/DiffGui"

echo "=== DiffGUI Setup for Expanse ==="
echo "Work directory: $WORK_DIR"

# Step 1: Clone the repository
if [ -d "$DIFFGUI_DIR" ]; then
    echo "DiffGui directory exists, pulling latest..."
    cd "$DIFFGUI_DIR" && git pull
else
    echo "Cloning DiffGui repository..."
    cd "$WORK_DIR"
    git clone https://github.com/QiaoyuHu89/DiffGui.git
fi

cd "$DIFFGUI_DIR"

# Step 2: Create conda environment
# The env.yml pins Python 3.7 + PyTorch 1.13 + CUDA 11.6
# On Expanse we need to adapt for the available CUDA modules

echo ""
echo "=== Creating conda environment ==="
echo "This will take 15-30 minutes..."

# Load modules
module purge
module load cpu/0.17.3b
module load gpu/0.17.3b
module load anaconda3/2021.05

# Create a trimmed env - the full env.yml has too many pinned system libs
# We'll install the key packages manually for compatibility with Expanse
conda create -n diffgui python=3.7 -y 2>/dev/null || true
conda activate diffgui

# Core dependencies
conda install pytorch=1.13.1 pytorch-cuda=11.6 -c pytorch -c nvidia -y
conda install pyg=2.3.1 -c pyg -y
conda install rdkit=2018.09.1 -c conda-forge -y
conda install openbabel=3.1.0 -c conda-forge -y
conda install biopython=1.78 scipy=1.7.3 scikit-learn=1.0.2 -y
conda install pandas=1.2.5 matplotlib seaborn networkx pyyaml easydict tqdm -y
conda install h5py python-lmdb -y

# Pip dependencies
pip install meeko==0.1.dev3 pdb2pqr vina==1.2.2
python -m pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3
pip install diffusers==0.21.4

# PyTorch extensions (torch_cluster, torch_scatter)
# These need to match PyTorch 1.13 + CUDA 11.6
pip install torch-scatter torch-cluster -f https://data.pyg.org/whl/torch-1.13.0+cu116.html

echo ""
echo "=== Downloading pretrained checkpoints ==="
echo "You need to manually download from Google Drive:"
echo "  https://drive.google.com/drive/folders/1pQk1FASCnCLjYRd7yc17WfctoHR50s2r"
echo ""
echo "Required files:"
echo "  - trained.pt          -> $DIFFGUI_DIR/ckpt/trained.pt"
echo "  - bond_trained.pt     -> $DIFFGUI_DIR/ckpt/bond_trained.pt"
echo ""
echo "On your local machine, use gdown or browser to download, then scp to Expanse:"
echo "  scp trained.pt bond_trained.pt $USER@login.expanse.sdsc.edu:$DIFFGUI_DIR/ckpt/"

mkdir -p "$DIFFGUI_DIR/ckpt"
mkdir -p "$DIFFGUI_DIR/outputs"

echo ""
echo "=== Setup Complete ==="
echo ""
echo "To use DiffGUI:"
echo "  1. Download checkpoints (see above)"
echo "  2. Prepare your pocket PDB (e.g., GID1 binding pocket)"
echo "  3. Edit configs/sample/sample.yml to point to your pocket"
echo "  4. Submit: sbatch submit_diffgui_sample.sh"
echo ""
echo "For scaffold-constrained generation:"
echo "  - Set gen_mode: frag_cond (keeps scaffold fixed, generates around it)"
echo "  - Provide scaffold as ligand: path/to/scaffold.sdf"
echo "  - Set frag: path/to/fragment_indices.txt"
