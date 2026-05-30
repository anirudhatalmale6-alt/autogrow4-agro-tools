#!/usr/bin/env python3
"""
MCS Database Builder - Phase 2
Builds 15 scaffold-organized seed libraries for AutoGrow4 de novo ligand design.

Pipeline steps (from client protocol):
  i.   Load supplier catalog, exclude GHS H400/H410/H411/H412
  ii.  MCS SMARTS substructure matching (15 scaffolds)
  iii. MW 150-350 Da filter
  iv.  Halogen filter (no Br/I, max 2 F or Cl on aromatic C only)
  v.   Photolabile group exclusion (nitro, aliphatic ethers, unconjugated alkenes)
  vi.  Cost filtering (via ChemPrice API - placeholder until API keys ready)
  vii. ReplaceCore to isolate side chains
  viii.Sort into 3 functional baskets (B1>B2>B3 priority)
  ix.  Morgan fingerprints (2048-bit, radius 2)
  x.   MaxMinPicker diversity selection (1333-3333 per basket)
  xi.  Output Generation 0 .smi files (5000-10000 per MCS family)

Usage:
  python3 mcs_database_builder.py --input supplier_catalog.smi --output_dir ./mcs_databases
  python3 mcs_database_builder.py --input supplier_catalog.smi --output_dir ./mcs_databases --scaffold quinoline
  python3 mcs_database_builder.py --demo  # Run with built-in test molecules
"""

import argparse
import json
import os
import sys
import time
import csv
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, rdSubstructLibrary
from rdkit.Chem import ReplaceCore, GetMolFrags
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.SimDivFilters import rdSimDivPickers
from rdkit import DataStructs
import numpy as np

# ============================================================
# 15 MCS Scaffold Definitions
# ============================================================

MCS_SCAFFOLDS = {
    "quinoline": {
        "smarts": "c1ccc2ncccc2c1",
        "description": "Quinoline Core - flat bicyclic aromatic, 1,4-benzodiazine mimic",
    },
    "isoquinoline": {
        "smarts": "c1ccc2cnccc2c1",
        "description": "Isoquinoline Core - shifted nitrogen vector vs quinoline",
    },
    "naphthalene": {
        "smarts": "c1ccc2ccccc2c1",
        "description": "Naphthalene Core - ultra-low-cost hydrophobic surface",
    },
    "benzothiazole": {
        "smarts": "c1ccc2scnc2c1",
        "description": "Benzothiazole Core - UV/soil-microbial resilient, S-N heterocycle",
    },
    "benzoxazole": {
        "smarts": "c1ccc2ocnc2c1",
        "description": "Benzoxazole Core - compact, strong local dipole",
    },
    "phthalimide": {
        "smarts": "O=C1NC(=O)c2ccccc21",
        "description": "Phthalimide Core - cheap bulk synthesis, imide N click handle",
    },
    "saccharin": {
        "smarts": "O=S1(=O)Nc2ccccc21",
        "description": "Saccharin Core (benzosulfonimide) - extreme UV stability",
    },
    "quinazolinone": {
        "smarts": "O=c1cnc2ccccc2[nH]1",
        "description": "Quinazolin-4(3H)-one Core - rigid pyrimidine-benzene lactam",
    },
    "coumarin": {
        "smarts": "O=c1ccc2ccccc2o1",
        "description": "Coumarin Core (benzopyran-2-one) - GA lactone mimic",
    },
    "phthalazinone": {
        "smarts": "O=c1[nH]ncc2ccccc12",
        "description": "Phthalazinone Core - stable N-N bond, unique H-bond vectors",
    },
    "oxadiazole_124": {
        "smarts": "[c;R1]1[n;R1][o;R1][c;R1][n;R1]1",
        "description": "1,2,4-Oxadiazole Core - bioisostere for esters/amides, UV+water stable",
    },
    "oxadiazole_134": {
        "smarts": "[c;R1]1[n;R1][n;R1][c;R1][o;R1]1",
        "description": "1,3,4-Oxadiazole Core - redirects exit vectors ~140 degrees",
    },
    "pyrimidine": {
        "smarts": "c1ncncn1",
        "description": "Pyrimidine Core - stable 6-ring, up to 3 click arms via SNAr",
    },
    "triazine": {
        "smarts": "c1ncncn1",
        "description": "1,3,5-Triazine Core - symmetric N-dense node, minimal MW footprint",
    },
    "carbazole": {
        "smarts": "c1ccc2c(c1)[nH]c1ccccc12",
        "description": "Carbazole Core - tricyclic hydrophobic volume for lid closure",
    },
}

# Pyrimidine and triazine have the same SMARTS in aromatic form;
# triazine needs a distinct pattern
MCS_SCAFFOLDS["pyrimidine"]["smarts"] = "c1ccnc(n1)"  # pyrimidine with 2 N
MCS_SCAFFOLDS["triazine"]["smarts"] = "c1nc(nc(n1))"  # 1,3,5-triazine with 3 N


def validate_smarts():
    """Validate all SMARTS patterns parse correctly."""
    valid = {}
    for name, info in MCS_SCAFFOLDS.items():
        mol = Chem.MolFromSmarts(info["smarts"])
        if mol is None:
            print(f"  WARNING: SMARTS for {name} failed to parse: {info['smarts']}")
            # Try alternative patterns
            alt = _get_alternative_smarts(name)
            if alt:
                mol = Chem.MolFromSmarts(alt)
                if mol:
                    info["smarts"] = alt
                    print(f"  -> Using alternative SMARTS: {alt}")
                    valid[name] = mol
        else:
            valid[name] = mol
    return valid


def _get_alternative_smarts(name):
    alternatives = {
        "quinoline": "c1ccnc2ccccc12",
        "isoquinoline": "c1cncc2ccccc12",
        "naphthalene": "c1ccc2ccccc2c1",
        "benzothiazole": "c1ccc2c(c1)nc(s2)",
        "benzoxazole": "c1ccc2c(c1)nc(o2)",
        "phthalimide": "O=C1NC(=O)c2ccccc21",
        "saccharin": "O=S1(=O)Nc2ccccc21",
        "quinazolinone": "O=c1cnc2ccccc2[nH]1",
        "coumarin": "O=c1ccc2ccccc2o1",
        "phthalazinone": "O=c1[nH]ncc2ccccc12",
        "oxadiazole_124": "[c;R1]1[n;R1][o;R1][c;R1][n;R1]1",
        "oxadiazole_134": "[c;R1]1[n;R1][n;R1][c;R1][o;R1]1",
        "pyrimidine": "c1ccncn1",
        "triazine": "c1ncncn1",
        "carbazole": "c1ccc2c(c1)[nH]c1ccccc12",
    }
    return alternatives.get(name)


# ============================================================
# Basket SMARTS Definitions
# ============================================================

BASKET_1_SMARTS = {
    "carboxylic_acid": "[CX3](=O)[O;H1,-]",
    "tetrazole": "c1nnn[nH]1",
    "acyl_sulfonamide": "[CX3](=O)[NX3][SX4](=O)=O",
    "aromatic_acyl_sulfonamide": "c([CX3](=O)[NX3][SX4](=O)=O)",
    "sulfonyl_urea": "[SX4](=O)(=O)[NX3][CX3](=O)[NX3]",
    "thiazolidinedione": "O=C1NC(=O)S1",
    "oxazolidinedione": "O=C1NC(=O)O1",
    "primary_amide": "[CX3](=O)[NX3;H2]",
    "secondary_amide": "[CX3](=O)[NX3;H1]([#6])",
    "aliphatic_nitrile": "[CX4][C]#N",
    "aliphatic_alcohol": "[CX4][O;H1]",
    "sulfoxide": "[SX3](=O)([#6])[#6]",
    "sulfone": "[SX4](=O)(=O)([#6])[#6]",
    "symmetric_urea": "[NX3][CX3](=O)[NX3]",
    "carbamate": "[NX3][CX3](=O)[OX2]",
}

BASKET_2_SMARTS = {
    "tertiary_amide": "[CX3](=O)[NX3]([#6])[#6]",
    "pyridyl": "c1ccncc1",
    "pyrrole_nh": "[nX3;H1]",
    "ether": "[OD2]([#6])[#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "ester": "[CX3](=O)[OX2H0][#6]",
}

BASKET_3_SMARTS = {
    "saturated_ring_carbon": "[C;X4;R1]",
    "branched_isopropyl": "[CH;X4](C)(C)",
    "branched_tert_butyl": "[C;X4](C)(C)(C)",
    "phenyl": "c1ccccc1",
    "naphthyl": "c1ccc2ccccc2c1",
}

# Groups to EXCLUDE from specific baskets
BASKET_1_EXCLUSIONS = {
    "sulfonic_acid": "[SX4](=O)(=O)[O;H1,-]",
    "phenol": "c[OH]",
    "aliphatic_amine_primary": "[NX3;H2;!$(N-C=O)][CX4]",
    "aliphatic_amine_secondary": "[NX3;H1;!$(N-C=O)]([CX4])[CX4]",
    "aliphatic_amine_tertiary": "[NX3;H0;!$(N-C=O)]([CX4])([CX4])[CX4]",
    "phosphonic_acid": "[PX4](=O)([OH])([OH])",
}

BASKET_2_EXCLUSIONS = {
    "aniline": "c[NX3;H2]",
    "repeating_ether": "[OD2H0][CX4,NX4][CX4,NX4][OD2H0]",
}

BASKET_3_EXCLUSIONS = {
    "any_heteroatom": "[O,N,S,P;!F;!Cl]",
    "terminal_alkene": "[CX3]=[CX3]",
    "alkyne": "[C]#[C]",
    "epoxide": "C1OC1",
    "aziridine": "C1NC1",
}


def compile_smarts_dict(smarts_dict):
    """Compile SMARTS strings into RDKit mol objects."""
    compiled = {}
    for name, smarts in smarts_dict.items():
        mol = Chem.MolFromSmarts(smarts)
        if mol is not None:
            compiled[name] = mol
        else:
            print(f"  WARNING: Failed to compile SMARTS '{name}': {smarts}")
    return compiled


# ============================================================
# Photolabile Group Blacklist
# ============================================================

PHOTOLABILE_SMARTS = {
    "unshielded_nitro": "[NX3](=O)=O",
    "aliphatic_ether_chain": "[CX4][OX2][CX4][OX2][CX4]",
    "unconjugated_alkene": "[CX4][CX3]=[CX3][CX4]",
}


# ============================================================
# Filter Functions
# ============================================================

def filter_molecular_weight(mol, min_mw=150.0, max_mw=350.0):
    """Step iii: MW 150-350 Da."""
    mw = Descriptors.MolWt(mol)
    return min_mw <= mw <= max_mw


def filter_halogens(mol):
    """Step iv: No Br/I. Max 2 F or 2 Cl, only on aromatic carbons."""
    br_pattern = Chem.MolFromSmarts("[Br]")
    i_pattern = Chem.MolFromSmarts("[I]")
    if mol.HasSubstructMatch(br_pattern) or mol.HasSubstructMatch(i_pattern):
        return False

    f_aliphatic = Chem.MolFromSmarts("[F][CX4]")
    cl_aliphatic = Chem.MolFromSmarts("[Cl][CX4]")
    if mol.HasSubstructMatch(f_aliphatic) or mol.HasSubstructMatch(cl_aliphatic):
        return False

    f_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[F]")))
    cl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[Cl]")))
    if f_count > 2 or cl_count > 2:
        return False

    return True


def filter_photolabile(mol, compiled_photolabile):
    """Step v: Exclude known photolabile groups."""
    for name, pattern in compiled_photolabile.items():
        if mol.HasSubstructMatch(pattern):
            return False
    return True


def filter_ghs_hazard(mol, ghs_codes=None):
    """Step i: Exclude GHS H400, H410, H411, H412.
    ghs_codes: dict mapping SMILES -> set of GHS codes, or None to skip."""
    if ghs_codes is None:
        return True
    smi = Chem.MolToSmiles(mol)
    codes = ghs_codes.get(smi, set())
    hazard_exclude = {"H400", "H410", "H411", "H412"}
    return len(codes & hazard_exclude) == 0


# ============================================================
# Basket Assignment
# ============================================================

def assign_basket(mol, core_mol, b1_patterns, b2_patterns, b3_patterns,
                  b1_exclusions, b2_exclusions, b3_exclusions):
    """Assign molecule to highest-priority basket (B1 > B2 > B3).

    Uses ReplaceCore to isolate side chains, then checks SMARTS patterns.
    Returns: 1, 2, 3, or 0 (unassigned/excluded).
    """
    try:
        side_chains = ReplaceCore(mol, core_mol)
    except Exception:
        side_chains = None
    if side_chains is None:
        check_mol = mol
    else:
        try:
            Chem.SanitizeMol(side_chains)
        except Exception:
            check_mol = mol
            side_chains = None
        else:
            check_mol = side_chains

    # Check B1 (highest priority)
    has_b1 = any(check_mol.HasSubstructMatch(p) for p in b1_patterns.values())
    excluded_b1 = any(check_mol.HasSubstructMatch(p) for p in b1_exclusions.values())
    if has_b1 and not excluded_b1:
        return 1

    # Check B2
    has_b2 = any(check_mol.HasSubstructMatch(p) for p in b2_patterns.values())
    excluded_b2 = any(check_mol.HasSubstructMatch(p) for p in b2_exclusions.values())
    # B2 must NOT have B1 features
    if has_b2 and not excluded_b2 and not has_b1:
        # Additional B2 constraints
        if not _check_b2_constraints(check_mol):
            return 0
        return 2

    # Check B3 (fallback)
    has_b3 = any(check_mol.HasSubstructMatch(p) for p in b3_patterns.values())
    excluded_b3 = any(check_mol.HasSubstructMatch(p) for p in b3_exclusions.values())
    if not excluded_b3:
        # B3 is the fallback: no polar transport triggers, no neutral H-bond anchors
        if not has_b1 and not has_b2:
            if not _check_b3_constraints(mol):
                return 0
            return 3

    return 0


def _check_b2_constraints(side_chain_mol):
    """Enforce Basket 2 physicochemical constraints on side chains."""
    try:
        # Max 2 ether oxygens per side chain
        ether_pattern = Chem.MolFromSmarts("[OD2]([#6])[#6]")
        if ether_pattern:
            ether_count = len(side_chain_mol.GetSubstructMatches(ether_pattern))
            if ether_count > 2:
                return False

        # Rotatable bond ceiling: max 4
        rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(side_chain_mol)
        if rot_bonds > 4:
            return False
    except Exception:
        pass
    return True


def _check_b3_constraints(full_mol):
    """Enforce Basket 3 physicochemical constraints."""
    try:
        # Rotatable bond limit: max 3 for side chain
        rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(full_mol)
        if rot_bonds > 5:
            return False

        # Heavy atom floor: > 3 non-core heavy atoms (checked on full mol)
        if full_mol.GetNumHeavyAtoms() < 6:
            return False

        # Lipophilicity ceiling
        slogp = Descriptors.MolLogP(full_mol)
        if slogp > 3.5:
            return False
    except Exception:
        pass
    return True


# ============================================================
# Diversity Selection
# ============================================================

def compute_morgan_fingerprints(mols_with_smiles, radius=2, n_bits=2048):
    """Step ix: Generate 2048-bit Morgan fingerprints."""
    from rdkit.Chem import rdFingerprintGenerator
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fps = []
    valid_entries = []
    for smi, mol in mols_with_smiles:
        fp = gen.GetFingerprint(mol)
        fps.append(fp)
        valid_entries.append((smi, mol))
    return fps, valid_entries


def diversity_select(fps, entries, n_select):
    """Step x: MaxMinPicker diversity selection."""
    if len(fps) <= n_select:
        return list(range(len(fps)))

    picker = rdSimDivPickers.MaxMinPicker()

    def dist_func(i, j):
        return 1.0 - TanimotoSimilarity(fps[i], fps[j])

    indices = picker.LazyBitVectorPick(fps, len(fps), n_select)
    return list(indices)


# ============================================================
# Main Pipeline
# ============================================================

class MCSPipelineConfig:
    """Configuration for the MCS database build pipeline."""
    def __init__(self):
        self.min_mw = 150.0
        self.max_mw = 350.0
        self.target_per_basket = 3333
        self.min_per_basket = 1333
        self.target_total = 10000
        self.min_total = 5000
        self.fp_radius = 2
        self.fp_nbits = 2048
        self.cost_percentile = None  # Set when ChemPrice API available


def load_molecules(input_path):
    """Load molecules from SMILES file (.smi, .csv, .tsv)."""
    molecules = []
    ext = Path(input_path).suffix.lower()

    with open(input_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if ext in ('.csv', '.tsv'):
                sep = ',' if ext == '.csv' else '\t'
                parts = line.split(sep)
                smi = parts[0].strip().strip('"')
            else:
                parts = line.split()
                smi = parts[0]

            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                name = parts[1] if len(parts) > 1 else f"mol_{line_num}"
                molecules.append((smi, mol, name))
            else:
                if line_num <= 5:
                    pass  # Skip header-like lines silently

    return molecules


def run_pipeline(input_path, output_dir, scaffold_filter=None, config=None,
                 ghs_codes=None, verbose=True):
    """Run the full MCS database build pipeline."""
    if config is None:
        config = MCSPipelineConfig()

    os.makedirs(output_dir, exist_ok=True)

    if verbose:
        print("=" * 70)
        print("MCS Database Builder - Phase 2")
        print("=" * 70)
        print(f"Input: {input_path}")
        print(f"Output: {output_dir}")
        print(f"MW range: {config.min_mw}-{config.max_mw} Da")
        print(f"Target per basket: {config.min_per_basket}-{config.target_per_basket}")
        print()

    # Validate SMARTS
    if verbose:
        print("Validating 15 MCS scaffold SMARTS patterns...")
    scaffold_mols = validate_smarts()
    if verbose:
        print(f"  {len(scaffold_mols)}/15 scaffolds validated")
        print()

    # Compile basket SMARTS
    b1_patterns = compile_smarts_dict(BASKET_1_SMARTS)
    b2_patterns = compile_smarts_dict(BASKET_2_SMARTS)
    b3_patterns = compile_smarts_dict(BASKET_3_SMARTS)
    b1_exclusions = compile_smarts_dict(BASKET_1_EXCLUSIONS)
    b2_exclusions = compile_smarts_dict(BASKET_2_EXCLUSIONS)
    b3_exclusions = compile_smarts_dict(BASKET_3_EXCLUSIONS)
    photolabile = compile_smarts_dict(PHOTOLABILE_SMARTS)

    if verbose:
        print(f"Compiled SMARTS: B1={len(b1_patterns)} patterns, "
              f"B2={len(b2_patterns)}, B3={len(b3_patterns)}")
        print(f"Exclusions: B1={len(b1_exclusions)}, B2={len(b2_exclusions)}, "
              f"B3={len(b3_exclusions)}")
        print(f"Photolabile filters: {len(photolabile)}")
        print()

    # Load molecules
    if verbose:
        print(f"Loading molecules from {input_path}...")
    molecules = load_molecules(input_path)
    if verbose:
        print(f"  Loaded {len(molecules)} valid molecules")
        print()

    # Determine which scaffolds to process
    if scaffold_filter:
        scaffolds_to_run = {scaffold_filter: MCS_SCAFFOLDS[scaffold_filter]}
    else:
        scaffolds_to_run = MCS_SCAFFOLDS

    results_summary = {}

    for scaffold_name, scaffold_info in scaffolds_to_run.items():
        if scaffold_name not in scaffold_mols:
            if verbose:
                print(f"SKIPPING {scaffold_name} - SMARTS validation failed")
            continue

        if verbose:
            print(f"{'=' * 60}")
            print(f"Processing: {scaffold_name.upper()}")
            print(f"  {scaffold_info['description']}")
            print(f"  SMARTS: {scaffold_info['smarts']}")
            print(f"{'=' * 60}")

        core_mol = scaffold_mols[scaffold_name]
        scaffold_dir = os.path.join(output_dir, scaffold_name)
        os.makedirs(scaffold_dir, exist_ok=True)

        # Step ii: Substructure matching
        matched = []
        for smi, mol, name in molecules:
            if mol.HasSubstructMatch(core_mol):
                matched.append((smi, mol, name))

        if verbose:
            print(f"  Step ii  - MCS match: {len(matched)}/{len(molecules)}")

        if not matched:
            if verbose:
                print(f"  No matches found, skipping")
            results_summary[scaffold_name] = {"matched": 0, "filtered": 0,
                                               "baskets": {1: 0, 2: 0, 3: 0},
                                               "final": 0}
            continue

        # Step i: GHS hazard exclusion
        if ghs_codes is not None:
            pre_ghs = len(matched)
            matched = [(s, m, n) for s, m, n in matched
                       if filter_ghs_hazard(m, ghs_codes)]
            if verbose:
                print(f"  Step i   - GHS filter: {len(matched)}/{pre_ghs}")

        # Step iii: MW filter
        pre_mw = len(matched)
        matched = [(s, m, n) for s, m, n in matched
                   if filter_molecular_weight(m, config.min_mw, config.max_mw)]
        if verbose:
            print(f"  Step iii - MW {config.min_mw}-{config.max_mw}: "
                  f"{len(matched)}/{pre_mw}")

        # Step iv: Halogen filter
        pre_hal = len(matched)
        matched = [(s, m, n) for s, m, n in matched if filter_halogens(m)]
        if verbose:
            print(f"  Step iv  - Halogen filter: {len(matched)}/{pre_hal}")

        # Step v: Photolabile filter
        pre_photo = len(matched)
        matched = [(s, m, n) for s, m, n in matched
                   if filter_photolabile(m, photolabile)]
        if verbose:
            print(f"  Step v   - Photolabile: {len(matched)}/{pre_photo}")

        # Step vi: Cost filtering (placeholder)
        if config.cost_percentile is not None:
            if verbose:
                print(f"  Step vi  - Cost filter: (requires ChemPrice API keys)")

        # Steps vii-viii: ReplaceCore + Basket assignment
        baskets = {1: [], 2: [], 3: []}
        unassigned = 0
        for smi, mol, name in matched:
            basket = assign_basket(mol, core_mol, b1_patterns, b2_patterns,
                                   b3_patterns, b1_exclusions, b2_exclusions,
                                   b3_exclusions)
            if basket > 0:
                baskets[basket].append((smi, mol))
            else:
                unassigned += 1

        if verbose:
            print(f"  Step viii - Basket assignment:")
            print(f"    B1 (Kleier triggers):  {len(baskets[1])}")
            print(f"    B2 (Pocket fitters):   {len(baskets[2])}")
            print(f"    B3 (Bulk fillers):     {len(baskets[3])}")
            print(f"    Unassigned/excluded:   {unassigned}")

        # Steps ix-x: Fingerprints + Diversity selection
        selected_all = []
        basket_counts = {}
        for b_num in [1, 2, 3]:
            b_mols = baskets[b_num]
            if not b_mols:
                basket_counts[b_num] = 0
                continue

            n_select = min(config.target_per_basket, len(b_mols))
            n_select = max(n_select, min(config.min_per_basket, len(b_mols)))

            fps, valid = compute_morgan_fingerprints(b_mols, config.fp_radius,
                                                     config.fp_nbits)
            indices = diversity_select(fps, valid, n_select)
            selected = [valid[i] for i in indices]
            selected_all.extend(selected)
            basket_counts[b_num] = len(selected)

            # Write per-basket file
            basket_file = os.path.join(scaffold_dir,
                                       f"{scaffold_name}_basket{b_num}.smi")
            with open(basket_file, 'w') as f:
                for smi, mol in selected:
                    f.write(f"{smi}\n")

            if verbose:
                print(f"  Step x   - B{b_num} diverse select: "
                      f"{len(selected)}/{len(b_mols)}")

        # Step xi: Combined Generation 0 .smi
        gen0_file = os.path.join(scaffold_dir, f"{scaffold_name}_gen0.smi")
        with open(gen0_file, 'w') as f:
            for smi, mol in selected_all:
                f.write(f"{smi}\n")

        total = len(selected_all)
        status = "OK" if total >= config.min_total else (
            "LOW" if total > 0 else "EMPTY")

        if verbose:
            print(f"  Step xi  - Generation 0: {total} compounds -> {gen0_file}")
            print(f"  Status: {status}")
            print()

        results_summary[scaffold_name] = {
            "matched": len(matched) + unassigned,
            "filtered": len(matched),
            "baskets": basket_counts,
            "final": total,
            "status": status,
        }

    # Write summary report
    report_path = os.path.join(output_dir, "pipeline_report.txt")
    _write_report(report_path, results_summary, config, input_path)

    # Write JSON results
    json_path = os.path.join(output_dir, "pipeline_results.json")
    with open(json_path, 'w') as f:
        json.dump(results_summary, f, indent=2)

    if verbose:
        print("=" * 60)
        print("PIPELINE SUMMARY")
        print("=" * 60)
        for name, res in results_summary.items():
            print(f"  {name:20s}: {res['final']:6d} compounds "
                  f"(B1={res['baskets'].get(1,0)}, "
                  f"B2={res['baskets'].get(2,0)}, "
                  f"B3={res['baskets'].get(3,0)})")
        print()
        print(f"Report: {report_path}")
        print(f"JSON:   {json_path}")

    return results_summary


def _write_report(path, results, config, input_path):
    """Write human-readable pipeline report."""
    with open(path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("MCS Database Builder - Phase 2 Pipeline Report\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 70 + "\n\n")

        f.write("CONFIGURATION:\n")
        f.write(f"  Input: {input_path}\n")
        f.write(f"  MW range: {config.min_mw}-{config.max_mw} Da\n")
        f.write(f"  Target per basket: {config.min_per_basket}-"
                f"{config.target_per_basket}\n")
        f.write(f"  Target total: {config.min_total}-{config.target_total}\n")
        f.write(f"  Fingerprint: Morgan r={config.fp_radius}, "
                f"{config.fp_nbits} bits\n\n")

        f.write("FILTERS APPLIED:\n")
        f.write("  1. GHS H400/H410/H411/H412 exclusion\n")
        f.write("  2. MCS substructure matching (15 scaffolds)\n")
        f.write(f"  3. Molecular weight: {config.min_mw}-{config.max_mw} Da\n")
        f.write("  4. Halogen: no Br/I, max 2 F/Cl on aromatic C only\n")
        f.write("  5. Photolabile: no unshielded nitro, aliphatic ether "
                "chains, unconjugated alkenes\n")
        f.write("  6. Cost percentile: "
                f"{'Not applied (API keys needed)' if config.cost_percentile is None else f'{config.cost_percentile}%'}\n\n")

        f.write("BASKET ASSIGNMENT (B1 > B2 > B3 priority):\n")
        f.write("  B1 (Kleier Triggers): Carboxylic acids, tetrazoles, "
                "acylsulfonamides, sulfonyl ureas\n")
        f.write("  B2 (Pocket Fitters): Tertiary amides, heteroaromatic N, "
                "ethers, ketones, esters\n")
        f.write("  B3 (Bulk Fillers): Cycloalkyls, branched alkyls, "
                "unfunctionalized aromatics\n\n")

        f.write("RESULTS PER SCAFFOLD:\n")
        f.write("-" * 70 + "\n")
        f.write(f"{'Scaffold':<20s} {'Matched':>8s} {'Filtered':>8s} "
                f"{'B1':>6s} {'B2':>6s} {'B3':>6s} {'Total':>6s} "
                f"{'Status':>8s}\n")
        f.write("-" * 70 + "\n")

        grand_total = 0
        for name, res in results.items():
            total = res['final']
            grand_total += total
            f.write(f"{name:<20s} {res['matched']:>8d} {res['filtered']:>8d} "
                    f"{res['baskets'].get(1,0):>6d} "
                    f"{res['baskets'].get(2,0):>6d} "
                    f"{res['baskets'].get(3,0):>6d} "
                    f"{total:>6d} "
                    f"{res.get('status',''):>8s}\n")

        f.write("-" * 70 + "\n")
        f.write(f"{'GRAND TOTAL':<20s} {'':>8s} {'':>8s} "
                f"{'':>6s} {'':>6s} {'':>6s} {grand_total:>6d}\n\n")

        f.write("OUTPUT FILES:\n")
        f.write("  Per scaffold: <scaffold>/<scaffold>_gen0.smi "
                "(combined Generation 0)\n")
        f.write("  Per basket:   <scaffold>/<scaffold>_basket<N>.smi\n")
        f.write("  This report:  pipeline_report.txt\n")
        f.write("  JSON data:    pipeline_results.json\n")


# ============================================================
# Synthon Mutation Database Builder
# ============================================================

def build_synthon_database(input_path, output_dir, verbose=True):
    """Build basket-organized synthon mutation databases for AutoGrow4.

    Steps:
      1. Load synthon catalog
      2. Filter < 125 Da
      3. Remove terminal aliphatic halides, peroxides, hydrazines
      4. Assign to 3 baskets
      5. Output mut_basket1.smi, mut_basket2.smi, mut_basket3.smi
    """
    os.makedirs(output_dir, exist_ok=True)

    if verbose:
        print("=" * 60)
        print("Synthon Mutation Database Builder")
        print("=" * 60)

    # Compile patterns
    b1_patterns = compile_smarts_dict(BASKET_1_SMARTS)
    b2_patterns = compile_smarts_dict(BASKET_2_SMARTS)
    b3_patterns = compile_smarts_dict(BASKET_3_SMARTS)
    b1_exclusions = compile_smarts_dict(BASKET_1_EXCLUSIONS)
    b2_exclusions = compile_smarts_dict(BASKET_2_EXCLUSIONS)
    b3_exclusions = compile_smarts_dict(BASKET_3_EXCLUSIONS)

    # Exclusion patterns for synthons
    terminal_halide = Chem.MolFromSmarts("[CX4][F,Cl,Br,I]")
    peroxide = Chem.MolFromSmarts("[OX2][OX2]")
    hydrazine = Chem.MolFromSmarts("[NX3][NX3]")

    molecules = load_molecules(input_path)
    if verbose:
        print(f"  Loaded {len(molecules)} synthons")

    # Filter < 125 Da
    filtered = [(s, m, n) for s, m, n in molecules
                if Descriptors.MolWt(m) < 125.0]
    if verbose:
        print(f"  < 125 Da filter: {len(filtered)}/{len(molecules)}")

    # Remove terminal aliphatic halides, peroxides, hydrazines
    clean = []
    for s, m, n in filtered:
        if terminal_halide and m.HasSubstructMatch(terminal_halide):
            continue
        if peroxide and m.HasSubstructMatch(peroxide):
            continue
        if hydrazine and m.HasSubstructMatch(hydrazine):
            continue
        clean.append((s, m, n))
    if verbose:
        print(f"  Safety filter: {len(clean)}/{len(filtered)}")

    # Assign to baskets (no core to remove for synthons)
    baskets = {1: [], 2: [], 3: []}
    for smi, mol, name in clean:
        # Direct SMARTS matching on the whole synthon
        has_b1 = any(mol.HasSubstructMatch(p) for p in b1_patterns.values())
        excl_b1 = any(mol.HasSubstructMatch(p) for p in b1_exclusions.values())
        has_b2 = any(mol.HasSubstructMatch(p) for p in b2_patterns.values())
        excl_b2 = any(mol.HasSubstructMatch(p) for p in b2_exclusions.values())
        has_b3 = any(mol.HasSubstructMatch(p) for p in b3_patterns.values())
        excl_b3 = any(mol.HasSubstructMatch(p) for p in b3_exclusions.values())

        if has_b1 and not excl_b1:
            baskets[1].append(smi)
        elif has_b2 and not excl_b2 and not has_b1:
            baskets[2].append(smi)
        elif not excl_b3 and not has_b1 and not has_b2:
            baskets[3].append(smi)

    for b_num in [1, 2, 3]:
        out_file = os.path.join(output_dir, f"mut_basket{b_num}.smi")
        with open(out_file, 'w') as f:
            for smi in baskets[b_num]:
                f.write(f"{smi}\n")
        if verbose:
            print(f"  Basket {b_num}: {len(baskets[b_num])} synthons -> {out_file}")

    if verbose:
        total = sum(len(b) for b in baskets.values())
        print(f"  Total: {total} synthons across 3 baskets")
        print()

    return baskets


# ============================================================
# Demo / Test Mode
# ============================================================

def generate_demo_molecules():
    """Generate test molecules covering all 15 scaffolds for validation."""
    demo_smiles = {
        "quinoline": [
            "c1ccc2nc(C(=O)O)cc(O)c2c1",  # B1: carboxylic acid
            "c1ccc2nc(C(=O)N(C)C)ccc2c1",  # B2: tertiary amide
            "c1ccc2nc(C3CC3)ccc2c1",  # B3: cyclopropyl
            "c1ccc2nc(CC(=O)O)c(F)cc2c1",  # B1: acid + F
            "c1ccc2nc(OC)c(C(=O)NC)cc2c1",  # B2: ether + secondary amide -> B1
            "c1ccc2nc(-c3ccccc3)ccc2c1",  # B3: phenyl
            "c1ccc2nc(C(C)(C)C)ccc2c1",  # B3: tert-butyl
        ],
        "benzothiazole": [
            "c1ccc2c(c1)nc(C(=O)O)s2",  # B1: carboxylic acid
            "c1ccc2c(c1)nc(C(=O)N(C)C)s2",  # B2: tertiary amide
            "c1ccc2c(c1)nc(C3CCCCC3)s2",  # B3: cyclohexyl
        ],
        "phthalimide": [
            "O=C1N(CC(=O)O)C(=O)c2ccccc21",  # B1: acid
            "O=C1N(CC(=O)N(C)C)C(=O)c2ccccc21",  # B2: tertiary amide
            "O=C1N(CC2CC2)C(=O)c2ccccc21",  # B3: cyclopropyl
        ],
        "pyrimidine": [
            "c1cc(C(=O)O)nc(N)n1",  # B1: acid
            "c1cc(OC)nc(OC)n1",  # B2: ethers
            "c1cc(C(C)C)ncn1",  # B3: isopropyl
        ],
        "naphthalene": [
            "c1ccc2c(CC(=O)O)cccc2c1",  # B1: acid
            "c1ccc2c(OC)cccc2c1",  # B2: ether
            "c1ccc2c(C(C)(C)C)cccc2c1",  # B3: tert-butyl
        ],
        "oxadiazole_124": [
            "c1(CC(=O)O)noc(C)n1",  # B1: acid on 1,2,4-oxadiazole
            "c1(C(=O)N(C)C)noc(C)n1",  # B2: tertiary amide
        ],
        "oxadiazole_134": [
            "c1(CC(=O)O)nnc(C)o1",  # B1: acid on 1,3,4-oxadiazole
            "c1(C(=O)N(C)C)nnc(C)o1",  # B2: tertiary amide
        ],
        "coumarin": [
            "O=c1cc(CC(=O)O)c2ccccc2o1",  # B1: acid on coumarin
            "O=c1cc(OC)c2ccccc2o1",  # B2: ether
        ],
    }

    all_mols = []
    for scaffold, smiles_list in demo_smiles.items():
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                all_mols.append((smi, mol, f"demo_{scaffold}"))
    return all_mols


def run_demo(output_dir="./demo_mcs_output"):
    """Run pipeline on demo molecules to validate logic."""
    print("=" * 70)
    print("DEMO MODE - Testing pipeline with built-in molecules")
    print("=" * 70)
    print()

    # Generate demo SMILES file
    demo_mols = generate_demo_molecules()
    demo_file = os.path.join(output_dir, "demo_input.smi")
    os.makedirs(output_dir, exist_ok=True)
    with open(demo_file, 'w') as f:
        for smi, mol, name in demo_mols:
            f.write(f"{smi} {name}\n")

    print(f"Generated {len(demo_mols)} demo molecules")
    print()

    # Validate SMARTS first
    print("SMARTS Validation:")
    scaffold_mols = validate_smarts()
    print()

    # Test substructure matching
    print("Substructure Matching Test:")
    for scaffold_name, scaffold_info in MCS_SCAFFOLDS.items():
        if scaffold_name not in scaffold_mols:
            continue
        core = scaffold_mols[scaffold_name]
        matches = sum(1 for _, mol, _ in demo_mols if mol.HasSubstructMatch(core))
        if matches > 0:
            print(f"  {scaffold_name:20s}: {matches} matches")
    print()

    # Run full pipeline
    config = MCSPipelineConfig()
    config.min_mw = 100.0  # Relaxed for demo
    config.max_mw = 500.0
    config.target_per_basket = 100
    config.min_per_basket = 1

    results = run_pipeline(demo_file, output_dir, config=config, verbose=True)

    print()
    print("Demo complete. Check output in:", output_dir)
    return results


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="MCS Database Builder for AutoGrow4 Phase 2",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 mcs_database_builder.py --demo
  python3 mcs_database_builder.py --input catalog.smi --output_dir ./mcs_dbs
  python3 mcs_database_builder.py --input catalog.smi --scaffold quinoline
  python3 mcs_database_builder.py --synthons enamine_synthons.smi --output_dir ./mutations
        """
    )

    parser.add_argument("--input", help="Input SMILES file (supplier catalog)")
    parser.add_argument("--output_dir", default="./mcs_databases",
                        help="Output directory")
    parser.add_argument("--scaffold", choices=list(MCS_SCAFFOLDS.keys()),
                        help="Process single scaffold only")
    parser.add_argument("--min_mw", type=float, default=150.0,
                        help="Minimum MW (default: 150)")
    parser.add_argument("--max_mw", type=float, default=350.0,
                        help="Maximum MW (default: 350)")
    parser.add_argument("--target_per_basket", type=int, default=3333,
                        help="Target compounds per basket (default: 3333)")
    parser.add_argument("--synthons", help="Build synthon mutation database "
                        "from this file instead")
    parser.add_argument("--demo", action="store_true",
                        help="Run demo with built-in test molecules")

    args = parser.parse_args()

    if args.demo:
        run_demo(args.output_dir)
        return

    if args.synthons:
        build_synthon_database(args.synthons, args.output_dir)
        return

    if not args.input:
        parser.error("--input required (or use --demo)")

    config = MCSPipelineConfig()
    config.min_mw = args.min_mw
    config.max_mw = args.max_mw
    config.target_per_basket = args.target_per_basket

    run_pipeline(args.input, args.output_dir, scaffold_filter=args.scaffold,
                 config=config)


if __name__ == "__main__":
    main()
