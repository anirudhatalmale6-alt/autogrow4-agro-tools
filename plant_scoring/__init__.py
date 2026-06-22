"""
Plant Systemic Bioavailability Scoring Pipeline
================================================
Five-gate multiplicative cascade predicting phloem mobility from SMILES.

Gates:
  1. Cuticle penetration (logKow/pp-LFER)
  2. Phloem loading (passive + carrier)
  3. Vascular transport (2-D Kleier surface)
  4. Unloading (pH-dependent + vacuolar)
  5. Metabolic stability (reactive-site alerts)

Architecture: Master Checkpoint v2 (Blocks 1A-2E)
Implementation: Phases A-J per Developer Implementation Plan

Requires: rdkit, pandas, openpyxl, numpy, requests
"""

__version__ = "0.1.0"

from plant_scoring.standardize import standardize_smiles, batch_standardize, fetch_identifiers
from plant_scoring.descriptors import compute_descriptors, detect_labile_groups, estimate_pka_groups
from plant_scoring.pipeline import process_compound, process_batch, run_demo
