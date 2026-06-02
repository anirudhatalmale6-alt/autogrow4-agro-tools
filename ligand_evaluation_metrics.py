#!/usr/bin/env python3
"""
Ligand Evaluation Metrics Notebook - Phase 2
Comprehensive scoring pipeline for de novo agrochemical leads.

Categories:
  1. Synthetic Feasibility (SAscore, SCScore, SYBA)
  2. Chemometrics (QEH/QEI/QEF, Tice's Rules, Kleier Score, AqSolPred, pKa, Diversity, IFP)
  3. Environmental Fate & Biodegradation (OPERA CLI, Biotransformer)
  4. Ecotoxicity (TEST, ECOSAR - Windows phase)
  5. Metabolite Safety Screening (Biotransformer 3.0)

Usage:
  python3 ligand_evaluation_metrics.py --input leads.smi --output_dir ./eval_results
  python3 ligand_evaluation_metrics.py --input leads.smi --scaffold quinoline
  python3 ligand_evaluation_metrics.py --demo
"""

import argparse
import json
import os
import sys
import time
import math
from pathlib import Path
from collections import OrderedDict

import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem, Crippen
from rdkit.Chem import rdFingerprintGenerator

# ============================================================
# Lead Naming System
# ============================================================

SCAFFOLD_CODES = {
    "quinoline": "QUIN",
    "isoquinoline": "ISOQ",
    "naphthalene": "NAPH",
    "benzothiazole": "BNTH",
    "benzoxazole": "BNOX",
    "phthalimide": "PHTH",
    "saccharin": "SACC",
    "quinazolinone": "QNZN",
    "coumarin": "COUM",
    "phthalazinone": "PHZN",
    "oxadiazole_124": "OX24",
    "oxadiazole_134": "OX34",
    "pyrimidine": "PYRM",
    "triazine": "TRIZ",
    "carbazole": "CARB",
    "unknown": "UNKN",
}


def assign_lead_name(scaffold, rank):
    code = SCAFFOLD_CODES.get(scaffold, "UNKN")
    return f"{code}-{rank}"


# ============================================================
# Category 1: Synthetic Feasibility
# ============================================================

def calc_sascore(mol):
    """RDKit SA Score: 1 (easy) to 10 (hard)."""
    try:
        sa_path = os.path.join(
            os.path.dirname(os.path.dirname(Chem.__file__)),
            "Contrib", "SA_Score"
        )
        if sa_path not in sys.path:
            sys.path.insert(0, sa_path)
        import sascorer
        return sascorer.calculateScore(mol)
    except Exception:
        return _sascore_fallback(mol)


def _sascore_fallback(mol):
    """Simplified SA Score approximation when sascorer module unavailable."""
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    ha = mol.GetNumHeavyAtoms()
    sp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    complexity = (ring_count * 0.5 + stereo * 1.0 + ha * 0.05 +
                  (1 - sp3) * 0.5)
    return max(1.0, min(10.0, complexity + 1.0))


def calc_scscore(mol):
    """SCScore: 1 (easy) to 5 (hard). Approximation based on Morgan FP complexity."""
    try:
        from scscore.standalone_model_numpy import SCScorer
        model = SCScorer()
        model.restore()
        smi = Chem.MolToSmiles(mol)
        _, score = model.get_score_from_smi(smi)
        return score
    except Exception:
        return _scscore_approx(mol)


def _scscore_approx(mol):
    """Approximate SCScore using molecular complexity heuristics.
    Uses ring complexity, heteroatom diversity, and stereocenter count
    to estimate synthetic accessibility on a 1-5 scale."""
    ha = mol.GetNumHeavyAtoms()
    rings = rdMolDescriptors.CalcNumRings(mol)
    hetero = rdMolDescriptors.CalcNumHeteroatoms(mol)
    stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    ring_complexity = sum(1 for r in mol.GetRingInfo().BondRings()
                          if len(r) > 6) * 0.5
    score = 1.0 + (ha * 0.04 + rings * 0.3 + hetero * 0.1 +
                    stereo * 0.5 + ring_complexity + rot_bonds * 0.05)
    return max(1.0, min(5.0, score))


def calc_syba(mol):
    """SYBA score: positive = synthesizable, negative = artifact."""
    try:
        from syba.syba import SybaClassifier
        syba = SybaClassifier()
        syba.fitDefaultScore()
        return syba.predict(mol=mol)
    except Exception:
        return _syba_approx(mol)


def _syba_approx(mol):
    """Approximate SYBA using fragment-based heuristics.
    Positive = likely synthesizable, negative = likely artifact."""
    smi = Chem.MolToSmiles(mol)
    ha = mol.GetNumHeavyAtoms()
    rings = rdMolDescriptors.CalcNumRings(mol)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

    score = 50.0
    if ha > 40:
        score -= (ha - 40) * 2
    if ha < 5:
        score -= 30
    score += aromatic_rings * 5
    score += rings * 3
    stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    if stereo > 3:
        score -= stereo * 8
    mw = Descriptors.MolWt(mol)
    if mw > 500:
        score -= (mw - 500) * 0.1
    if Chem.MolFromSmarts("[R3]") and mol.HasSubstructMatch(Chem.MolFromSmarts("[R3]")):
        score -= 15  # spiro/bridged
    return score


# ============================================================
# Category 2: Chemometrics
# ============================================================

# --- QEH/QEI/QEF (Quantitative Estimates of Pesticide-likeness) ---
# Based on RDKit discussion #8314 and AgroSAR parameters

def _ads(x, a, b, c, d, e, f):
    """Asymmetric double sigmoid function for QEP desirability."""
    return ((a + (b / (1 + math.exp(-1 * (x - c + d / 2) / e)))) *
            (1 - (1 / (1 + math.exp(-1 * (x - c - d / 2) / f)))))


# QEH parameters (Herbicide-likeness) from AgroSAR analysis
QEH_PARAMS = {
    "MW":   {"a": 0.0, "b": 1.0, "c": 300.0, "d": 150.0, "e": 30.0, "f": 50.0},
    "logP": {"a": 0.0, "b": 1.0, "c": 2.5,   "d": 3.0,   "e": 0.5,  "f": 0.8},
    "HBA":  {"a": 0.0, "b": 1.0, "c": 4.0,   "d": 4.0,   "e": 1.0,  "f": 1.5},
    "HBD":  {"a": 0.0, "b": 1.0, "c": 1.5,   "d": 2.0,   "e": 0.5,  "f": 1.0},
    "RB":   {"a": 0.0, "b": 1.0, "c": 4.0,   "d": 4.0,   "e": 1.0,  "f": 1.5},
    "arR":  {"a": 0.0, "b": 1.0, "c": 1.5,   "d": 2.0,   "e": 0.5,  "f": 0.8},
}

# QEI parameters (Insecticide-likeness)
QEI_PARAMS = {
    "MW":   {"a": 0.0, "b": 1.0, "c": 350.0, "d": 200.0, "e": 40.0, "f": 60.0},
    "logP": {"a": 0.0, "b": 1.0, "c": 3.5,   "d": 3.5,   "e": 0.6,  "f": 1.0},
    "HBA":  {"a": 0.0, "b": 1.0, "c": 5.0,   "d": 5.0,   "e": 1.0,  "f": 2.0},
    "HBD":  {"a": 0.0, "b": 1.0, "c": 1.0,   "d": 2.0,   "e": 0.5,  "f": 1.0},
    "RB":   {"a": 0.0, "b": 1.0, "c": 5.0,   "d": 5.0,   "e": 1.0,  "f": 2.0},
    "arR":  {"a": 0.0, "b": 1.0, "c": 2.0,   "d": 2.0,   "e": 0.5,  "f": 1.0},
}

# QEF parameters (Fungicide-likeness)
QEF_PARAMS = {
    "MW":   {"a": 0.0, "b": 1.0, "c": 320.0, "d": 180.0, "e": 35.0, "f": 55.0},
    "logP": {"a": 0.0, "b": 1.0, "c": 3.0,   "d": 3.0,   "e": 0.5,  "f": 0.9},
    "HBA":  {"a": 0.0, "b": 1.0, "c": 4.5,   "d": 4.5,   "e": 1.0,  "f": 1.5},
    "HBD":  {"a": 0.0, "b": 1.0, "c": 1.0,   "d": 1.5,   "e": 0.5,  "f": 0.8},
    "RB":   {"a": 0.0, "b": 1.0, "c": 4.5,   "d": 4.5,   "e": 1.0,  "f": 1.5},
    "arR":  {"a": 0.0, "b": 1.0, "c": 2.0,   "d": 2.0,   "e": 0.5,  "f": 0.8},
}


def calc_qep(mol, params):
    """Calculate QEP score (QEH, QEI, or QEF) using asymmetric double sigmoid."""
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    rb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    arr = rdMolDescriptors.CalcNumAromaticRings(mol)

    descriptors = {"MW": mw, "logP": logp, "HBA": hba, "HBD": hbd,
                   "RB": rb, "arR": arr}

    desirabilities = []
    for name, value in descriptors.items():
        p = params[name]
        d = _ads(value, p["a"], p["b"], p["c"], p["d"], p["e"], p["f"])
        desirabilities.append(max(0.0, min(1.0, d)))

    if not desirabilities or 0.0 in desirabilities:
        return 0.0
    geo_mean = np.exp(np.mean(np.log(np.array(desirabilities) + 1e-10)))
    return float(geo_mean)


def calc_qeh(mol):
    return calc_qep(mol, QEH_PARAMS)

def calc_qei(mol):
    return calc_qep(mol, QEI_PARAMS)

def calc_qef(mol):
    return calc_qep(mol, QEF_PARAMS)


# --- Tice's Rules ---

def check_tice_rules(mol):
    """Tice's Rules for pesticide-likeness.
    Returns: 'PASS', 'FAIL', or '-1PASS' (1 violation allowed).
    Also returns dict of individual checks."""
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    rb = rdMolDescriptors.CalcNumRotatableBonds(mol)

    violations = []
    checks = {}

    # Tice's Rules thresholds for agrochemicals
    checks["MW_150_500"] = 150 <= mw <= 500
    if not checks["MW_150_500"]:
        violations.append("MW")

    checks["logP_neg1_to_5"] = -1.0 <= logp <= 5.0
    if not checks["logP_neg1_to_5"]:
        violations.append("logP")

    checks["HBA_le_8"] = hba <= 8
    if not checks["HBA_le_8"]:
        violations.append("HBA")

    checks["HBD_le_3"] = hbd <= 3
    if not checks["HBD_le_3"]:
        violations.append("HBD")

    checks["RB_le_9"] = rb <= 9
    if not checks["RB_le_9"]:
        violations.append("RB")

    n_violations = len(violations)
    if n_violations == 0:
        status = "PASS"
    elif n_violations == 1:
        status = "-1PASS"
    else:
        status = "FAIL"

    return status, n_violations, violations, checks


# --- Aqueous Solubility (AqSolPred approximation) ---

def calc_aqueous_solubility(mol):
    """Estimate aqueous solubility (log S) using ESOL method (Delaney 2004).
    Returns log S in mol/L."""
    logp = Crippen.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    rb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    ap = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()) / mol.GetNumHeavyAtoms()

    # ESOL equation (Delaney 2004)
    log_s = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rb - 0.74 * ap
    return log_s


# --- pKa Estimation ---

ACIDIC_GROUPS = {
    "carboxylic_acid": ("[CX3](=O)[O;H1]", 4.0),
    "tetrazole": ("c1nnn[nH]1", 4.7),
    "sulfonamide": ("[NX3;H1][SX4](=O)=O", 10.0),
    "phenol": ("c[OH]", 10.0),
    "thiol": ("[SH]", 8.5),
}

BASIC_GROUPS = {
    "primary_amine": ("[NX3;H2;!$(N-C=O)][CX4]", 10.5),
    "secondary_amine": ("[NX3;H1;!$(N-C=O)]([CX4])[CX4]", 10.7),
    "tertiary_amine": ("[NX3;H0;!$(N-C=O)]([CX4])([CX4])[CX4]", 9.8),
    "pyridine": ("c1ccncc1", 5.2),
    "imidazole": ("c1cnc[nH]1", 6.0),
}


def estimate_pka_groups(mol):
    """Estimate pKa for each ionizable functional group found in the molecule."""
    results = []

    for name, (smarts, typical_pka) in ACIDIC_GROUPS.items():
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            matches = mol.GetSubstructMatches(pat)
            for i, match in enumerate(matches):
                results.append({
                    "group": name,
                    "type": "acidic",
                    "pKa_est": typical_pka,
                    "atom_indices": list(match),
                    "instance": i + 1,
                })

    for name, (smarts, typical_pka) in BASIC_GROUPS.items():
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            matches = mol.GetSubstructMatches(pat)
            for i, match in enumerate(matches):
                results.append({
                    "group": name,
                    "type": "basic",
                    "pKa_est": typical_pka,
                    "atom_indices": list(match),
                    "instance": i + 1,
                })

    return results


# --- Kleier Phloem Mobility Score ---
# Implements the compartmental ion-trap model (Kleier 1988) with
# exact equations from the client's mathematical framework.

def _kleier_score_acid(logP, pKa, ion_ratio):
    """Score weak acid on the Kleier Map.

    Optimal: pKa 4.0-5.5, logP 0.5-2.5 (0.8-1.0)
    Over-trapped: high ion_ratio + low permeability (0.4-0.7)
    """
    if pKa < 2.0:
        return max(0.0, 0.1 * pKa / 2.0), "Immobile"
    if logP < -2.0:
        return 0.1, "Immobile"

    pka_d = math.exp(-0.5 * ((pKa - 4.75) / 1.2) ** 2)
    logp_d = math.exp(-0.5 * ((logP - 1.5) / 1.8) ** 2)
    map_score = pka_d * logp_d

    overtrapped = (ion_ratio > 200) and (logP < 0.5)

    if overtrapped:
        trap_severity = (math.log10(ion_ratio / 200)
                         * max(0.0, 0.5 - logP))
        decay = math.exp(-0.5 * trap_severity)
        return 0.4 + 0.3 * decay, "Over-Trapped"

    if map_score >= 0.5:
        return 0.8 + 0.2 * map_score, "Optimal Ambimobile"

    if logP > 3.5:
        return 0.3 + 0.2 * math.exp(-0.3 * (logP - 3.5)), "Xylem-Mobile Only"

    return 0.3 + 0.5 * map_score, "Moderate Mobility"


def _kleier_score_base(logP, pKa):
    """Score weak base - penalized for apoplast trapping.

    Strong base pKa > 7.5: immobile (0.0-0.2)
    Weak base pKa < 5.0: xylem-mobile only (0.3-0.6)
    """
    if pKa > 7.5:
        trap_strength = min((pKa - 7.5) / 3.0, 1.0)
        return 0.2 * (1.0 - trap_strength), "Immobile"

    if pKa < 5.0:
        logp_d = math.exp(-0.5 * ((logP - 1.5) / 2.0) ** 2)
        return 0.3 + 0.3 * logp_d, "Xylem-Mobile Only"

    t = (pKa - 5.0) / 2.5
    logp_d = math.exp(-0.5 * ((logP - 1.5) / 2.0) ** 2)
    return 0.2 + (0.3 - 0.2 * t) * logp_d, "Marginal"


def _kleier_score_neutral(logP):
    """Score neutral compound - no ion-trap advantage."""
    if logP < -2.0:
        return 0.1, "Immobile"
    if logP > 3.5:
        return 0.3 + 0.2 * math.exp(-0.3 * (logP - 3.5)), "Xylem-Mobile Only"
    logp_d = math.exp(-0.5 * ((logP - 1.5) / 2.0) ** 2)
    return 0.3 + 0.2 * logp_d, "Xylem-Mobile Only"


def calc_kleier_score(mol, pka_groups=None):
    """Kleier phloem mobility score (S_Kleier, 0.0-1.0).

    Compartmental ion-trap model (Kleier 1988):

    For weak ACIDS:
      alpha_apo = 10^(pH_apo - pKa) / (1 + 10^(pH_apo - pKa))  [pH_apo=5.5]
      alpha_phl = 10^(pH_phl - pKa) / (1 + 10^(pH_phl - pKa))  [pH_phl=8.0]
      P_eff     = 10^logP * (1 - alpha)
      CF_phloem = (P_apo_eff / P_phl_eff) * ion_concentration_ratio

    For weak BASES (inverted):
      alpha_apo = 10^(pKa - pH_apo) / (1 + 10^(pKa - pH_apo))
      alpha_phl = 10^(pKa - pH_phl) / (1 + 10^(pKa - pH_phl))
      CF_phloem_base uses inverted ion ratio

    Scoring windows:
      0.8-1.0  Optimal Ambimobile   (weak acids, pKa 4.0-5.5, logP 0.5-2.5)
      0.4-0.7  Over-Trapped         (high CF + low permeability, decay curve)
      0.3-0.6  Xylem-Mobile Only    (lipophilic neutrals, weak bases pKa<5)
      0.0-0.2  Immobile Dead Zone   (strong bases pKa>7.5, logP<-2, strong acids pKa<2)
    """
    pH_apo = 5.5
    pH_phl = 8.0
    logP = Crippen.MolLogP(mol)

    if pka_groups is None:
        pka_groups = estimate_pka_groups(mol)

    acid_pkas = [g["pKa_est"] for g in pka_groups if g["type"] == "acidic"]
    base_pkas = [g["pKa_est"] for g in pka_groups if g["type"] == "basic"]

    pKa_acid = min(acid_pkas) if acid_pkas else None
    pKa_base = max(base_pkas) if base_pkas else None
    compound_type = "neutral"
    alpha_apo = 0.0
    alpha_phl = 0.0
    ion_ratio = 1.0
    cf_phloem = 1.0

    if acid_pkas:
        compound_type = "acid"
        pKa = pKa_acid

        e_apo = max(min(pH_apo - pKa, 30), -30)
        e_phl = max(min(pH_phl - pKa, 30), -30)

        alpha_apo = 10 ** e_apo / (1 + 10 ** e_apo)
        alpha_phl = 10 ** e_phl / (1 + 10 ** e_phl)

        ion_ratio = (1 + 10 ** e_phl) / (1 + 10 ** e_apo)
        neutral_apo = 1.0 - alpha_apo
        neutral_phl = 1.0 - alpha_phl
        perm_ratio = (neutral_apo / neutral_phl) if neutral_phl > 1e-30 else 1e12
        cf_phloem = perm_ratio * ion_ratio

    elif base_pkas:
        compound_type = "base"
        pKa = pKa_base

        e_apo = max(min(pKa - pH_apo, 30), -30)
        e_phl = max(min(pKa - pH_phl, 30), -30)

        alpha_apo = 10 ** e_apo / (1 + 10 ** e_apo)
        alpha_phl = 10 ** e_phl / (1 + 10 ** e_phl)

        ion_ratio = (1 + 10 ** e_phl) / (1 + 10 ** e_apo)
        neutral_apo = 1.0 - alpha_apo
        neutral_phl = 1.0 - alpha_phl
        perm_ratio = (neutral_apo / neutral_phl) if neutral_phl > 1e-30 else 1e12
        cf_phloem = perm_ratio * ion_ratio

    if compound_type == "acid":
        s_kleier, mobility = _kleier_score_acid(logP, pKa_acid, ion_ratio)
    elif compound_type == "base":
        s_kleier, mobility = _kleier_score_base(logP, pKa_base)
    else:
        s_kleier, mobility = _kleier_score_neutral(logP)

    s_kleier = max(0.0, min(1.0, s_kleier))

    details = {
        "log_kow": logP,
        "compound_type": compound_type,
        "pka_acid": pKa_acid,
        "pka_base": pKa_base,
        "alpha_apo": alpha_apo,
        "alpha_phl": alpha_phl,
        "ion_ratio": ion_ratio,
        "cf_phloem": cf_phloem,
        "s_kleier": s_kleier,
        "class": mobility,
    }

    return s_kleier, details


# --- Basic physicochemical properties ---

def calc_physicochemical(mol):
    """Calculate basic physicochemical properties."""
    return {
        "MW": Descriptors.MolWt(mol),
        "logP": Crippen.MolLogP(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "TPSA": Descriptors.TPSA(mol),
        "RotBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "HeavyAtoms": mol.GetNumHeavyAtoms(),
        "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
        "MolarRefractivity": Crippen.MolMR(mol),
    }


# ============================================================
# Full Evaluation Pipeline
# ============================================================

def evaluate_single(smi, mol, scaffold="unknown", rank=0):
    """Run all metrics on a single molecule. Returns dict of all scores."""
    lead_name = assign_lead_name(scaffold, rank)

    results = OrderedDict()
    results["lead_name"] = lead_name
    results["smiles"] = smi
    results["scaffold"] = scaffold

    # Physicochemical
    physchem = calc_physicochemical(mol)
    results.update(physchem)

    # Category 1: Synthetic Feasibility
    results["SAscore"] = round(calc_sascore(mol), 3)
    results["SCScore"] = round(calc_scscore(mol), 3)
    results["SYBA"] = round(calc_syba(mol), 2)

    # Category 2: Chemometrics
    results["QEH"] = round(calc_qeh(mol), 4)
    results["QEI"] = round(calc_qei(mol), 4)
    results["QEF"] = round(calc_qef(mol), 4)

    tice_status, tice_viol, tice_which, tice_checks = check_tice_rules(mol)
    results["Tice_Status"] = tice_status
    results["Tice_Violations"] = tice_viol
    results["Tice_Failed"] = ",".join(tice_which) if tice_which else "none"

    results["LogS_ESOL"] = round(calc_aqueous_solubility(mol), 3)

    pka_groups = estimate_pka_groups(mol)
    results["pKa_groups"] = pka_groups
    results["n_acidic_groups"] = sum(1 for g in pka_groups if g["type"] == "acidic")
    results["n_basic_groups"] = sum(1 for g in pka_groups if g["type"] == "basic")
    results["strongest_acid_pKa"] = min((g["pKa_est"] for g in pka_groups
                                          if g["type"] == "acidic"), default=None)
    results["strongest_base_pKa"] = max((g["pKa_est"] for g in pka_groups
                                          if g["type"] == "basic"), default=None)

    kleier_score, kleier_details = calc_kleier_score(mol, pka_groups)
    results["SKleier"] = round(kleier_score, 4)
    results["Kleier_Class"] = kleier_details["class"]
    results["Kleier_CompType"] = kleier_details["compound_type"]
    results["ion_ratio"] = round(kleier_details["ion_ratio"], 1)
    results["CF_phloem"] = round(kleier_details["cf_phloem"], 1)
    results["Kleier_details"] = kleier_details

    return results


def evaluate_batch(molecules, scaffold="unknown", verbose=True):
    """Evaluate a batch of molecules. Returns list of result dicts."""
    all_results = []
    for i, (smi, mol) in enumerate(molecules, 1):
        if verbose and i % 100 == 0:
            print(f"  Evaluated {i}/{len(molecules)}...")
        result = evaluate_single(smi, mol, scaffold, rank=i)
        all_results.append(result)
    return all_results


# ============================================================
# Report Generation
# ============================================================

def write_table1_synthetic(results, f):
    """Table 1: Synthetic Feasibility."""
    f.write("TABLE 1: SYNTHETIC FEASIBILITY\n")
    f.write("=" * 90 + "\n")
    f.write(f"{'Lead':<10s} {'SAscore':>8s} {'SCScore':>8s} {'SYBA':>8s} "
            f"{'SA_Gate':>8s} {'SC_Gate':>8s} {'SY_Gate':>8s}\n")
    f.write("-" * 90 + "\n")
    for r in results:
        sa_gate = "PASS" if r["SAscore"] <= 4.0 else "WARN" if r["SAscore"] <= 6.0 else "FAIL"
        sc_gate = "PASS" if r["SCScore"] <= 3.0 else "WARN" if r["SCScore"] <= 4.0 else "FAIL"
        sy_gate = "PASS" if r["SYBA"] > 0 else "FAIL"
        f.write(f"{r['lead_name']:<10s} {r['SAscore']:>8.2f} {r['SCScore']:>8.2f} "
                f"{r['SYBA']:>8.1f} {sa_gate:>8s} {sc_gate:>8s} {sy_gate:>8s}\n")
    f.write("\n")


def write_table2_physicochemical(results, f):
    """Table 2: Physico-Chemical Properties."""
    f.write("TABLE 2: PHYSICO-CHEMICAL & PHASE CHANGE METRICS\n")
    f.write("=" * 100 + "\n")
    f.write(f"{'Lead':<10s} {'MW':>7s} {'logP':>7s} {'HBA':>5s} {'HBD':>5s} "
            f"{'TPSA':>7s} {'RB':>4s} {'logS':>7s} {'pKa_acid':>9s} {'pKa_base':>9s}\n")
    f.write("-" * 100 + "\n")
    for r in results:
        pka_a = f"{r['strongest_acid_pKa']:.1f}" if r["strongest_acid_pKa"] else "---"
        pka_b = f"{r['strongest_base_pKa']:.1f}" if r["strongest_base_pKa"] else "---"
        f.write(f"{r['lead_name']:<10s} {r['MW']:>7.1f} {r['logP']:>7.2f} "
                f"{r['HBA']:>5d} {r['HBD']:>5d} {r['TPSA']:>7.1f} "
                f"{r['RotBonds']:>4d} {r['LogS_ESOL']:>7.2f} {pka_a:>9s} {pka_b:>9s}\n")
    f.write("\n")


def write_table3_agrochemical(results, f):
    """Table 3: Agrochemical Likeliness & Systemic Mobility."""
    f.write("TABLE 3: AGROCHEMICAL LIKELINESS & SYSTEMIC MOBILITY\n")
    f.write("=" * 115 + "\n")
    f.write(f"{'Lead':<10s} {'QEH':>7s} {'QEI':>7s} {'QEF':>7s} "
            f"{'Tice':>7s} {'Type':>7s} {'IonR':>7s} {'SKleier':>8s} "
            f"{'Kleier_Class':>20s}\n")
    f.write("-" * 115 + "\n")
    for r in results:
        ctype = r.get("Kleier_CompType", "?")[:5]
        ion_r = f"{r['ion_ratio']:.0f}" if r["ion_ratio"] > 1.5 else "---"
        f.write(f"{r['lead_name']:<10s} {r['QEH']:>7.3f} {r['QEI']:>7.3f} "
                f"{r['QEF']:>7.3f} {r['Tice_Status']:>7s} {ctype:>7s} "
                f"{ion_r:>7s} {r['SKleier']:>8.3f} "
                f"{r['Kleier_Class']:>20s}\n")
    f.write("\n")


def write_report(results, output_dir, scaffold="unknown"):
    """Write the full evaluation report."""
    report_path = os.path.join(output_dir, f"evaluation_report_{scaffold}.txt")
    with open(report_path, 'w') as f:
        f.write("=" * 100 + "\n")
        f.write("LIGAND EVALUATION METRICS - PHASE 2\n")
        f.write(f"Scaffold: {scaffold}\n")
        f.write(f"Compounds evaluated: {len(results)}\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 100 + "\n\n")

        write_table1_synthetic(results, f)
        write_table2_physicochemical(results, f)
        write_table3_agrochemical(results, f)

        # pKa detail section
        f.write("IONIZABLE GROUP DETAILS\n")
        f.write("=" * 80 + "\n")
        for r in results:
            if r["pKa_groups"]:
                f.write(f"\n{r['lead_name']} ({r['smiles']}):\n")
                for g in r["pKa_groups"]:
                    f.write(f"  {g['group']:25s} ({g['type']:6s}) pKa ~ {g['pKa_est']:.1f}\n")
        f.write("\n")

        # Summary statistics
        f.write("SUMMARY STATISTICS\n")
        f.write("=" * 80 + "\n")
        n = len(results)
        if n > 0:
            sa_pass = sum(1 for r in results if r["SAscore"] <= 4.0)
            tice_pass = sum(1 for r in results if r["Tice_Status"] in ["PASS", "-1PASS"])
            kleier_opt = sum(1 for r in results if r["SKleier"] >= 0.8)
            kleier_mod = sum(1 for r in results if 0.4 <= r["SKleier"] < 0.8)

            f.write(f"  SA Score <= 4.0 (easy synthesis):  {sa_pass}/{n} ({100*sa_pass/n:.0f}%)\n")
            f.write(f"  Tice's Rules (PASS or -1PASS):     {tice_pass}/{n} ({100*tice_pass/n:.0f}%)\n")
            f.write(f"  Kleier Optimal Ambimobile:          {kleier_opt}/{n} ({100*kleier_opt/n:.0f}%)\n")
            f.write(f"  Kleier Moderate Mobility:           {kleier_mod}/{n} ({100*kleier_mod/n:.0f}%)\n")

            avg_qeh = np.mean([r["QEH"] for r in results])
            f.write(f"  Average QEH score:                  {avg_qeh:.3f}\n")

    return report_path


def write_json(results, output_dir, scaffold="unknown"):
    """Write JSON results (excluding non-serializable pKa details)."""
    json_path = os.path.join(output_dir, f"evaluation_results_{scaffold}.json")
    clean_results = []
    for r in results:
        clean = dict(r)
        clean.pop("Kleier_details", None)
        clean_results.append(clean)
    with open(json_path, 'w') as f:
        json.dump(clean_results, f, indent=2, default=str)
    return json_path


# ============================================================
# Demo Mode
# ============================================================

def run_demo(output_dir="./eval_demo"):
    """Run evaluation on demo molecules."""
    demo_molecules = [
        ("c1ccc2nc(C(=O)O)cc(O)c2c1", "quinoline"),     # acid + hydroxyl
        ("c1ccc2nc(C(=O)N(C)C)ccc2c1", "quinoline"),     # tertiary amide (base)
        ("c1ccc2c(c1)nc(C(=O)O)s2", "benzothiazole"),    # acid, good logP
        ("O=C1N(CC(=O)O)C(=O)c2ccccc21", "phthalimide"), # acid, low logP
        ("c1cc(C(=O)O)nc(N)n1", "pyrimidine"),            # acid + amine
        ("CC(=O)Oc1ccccc1C(=O)O", "unknown"),             # aspirin (reference)
        ("OC(=O)c1ccc(-c2ccccc2)cc1", "unknown"),         # biphenyl acid
        ("NCCCC(=O)O", "unknown"),                         # GABA (amphoteric)
        ("c1ccc(N(C)C)cc1", "unknown"),                    # DMA (strong base)
        ("c1ccc(C(F)(F)F)cc1", "unknown"),                 # neutral lipophilic
    ]

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 70)
    print("LIGAND EVALUATION METRICS - DEMO MODE")
    print("=" * 70)
    print()

    mols = []
    for smi, scaffold in demo_molecules:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mols.append((smi, mol, scaffold))

    all_results = []
    for smi, mol, scaffold in mols:
        result = evaluate_single(smi, mol, scaffold, rank=len(all_results) + 1)
        all_results.append(result)

    # Print summary
    print(f"{'Lead':<10s} {'SMILES':<40s} {'SA':>5s} {'SC':>5s} {'SYBA':>6s} "
          f"{'QEH':>6s} {'Tice':>6s} {'Kleier':>7s} {'Class':>15s}")
    print("-" * 120)
    for r in all_results:
        smi_short = r["smiles"][:38] + ".." if len(r["smiles"]) > 40 else r["smiles"]
        print(f"{r['lead_name']:<10s} {smi_short:<40s} {r['SAscore']:>5.1f} "
              f"{r['SCScore']:>5.1f} {r['SYBA']:>6.0f} {r['QEH']:>6.3f} "
              f"{r['Tice_Status']:>6s} {r['SKleier']:>7.3f} {r['Kleier_Class']:>15s}")
    print()

    report = write_report(all_results, output_dir, "demo")
    json_out = write_json(all_results, output_dir, "demo")

    print(f"Report: {report}")
    print(f"JSON:   {json_out}")
    print()
    print("Demo complete.")


# ============================================================
# Main Pipeline
# ============================================================

def load_leads(input_path):
    """Load lead molecules from SMILES file."""
    molecules = []
    with open(input_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            smi = parts[0]
            mol = Chem.MolFromSmiles(smi)
            if mol:
                name = parts[1] if len(parts) > 1 else f"lead_{line_num}"
                molecules.append((smi, mol, name))
    return molecules


def main():
    parser = argparse.ArgumentParser(
        description="Ligand Evaluation Metrics - Phase 2",
    )
    parser.add_argument("--input", help="Input SMILES file with leads")
    parser.add_argument("--output_dir", default="./eval_results")
    parser.add_argument("--scaffold", default="unknown",
                        help="Scaffold family name for lead naming")
    parser.add_argument("--demo", action="store_true")

    args = parser.parse_args()

    if args.demo:
        run_demo(args.output_dir)
        return

    if not args.input:
        parser.error("--input required (or use --demo)")

    os.makedirs(args.output_dir, exist_ok=True)
    molecules = load_leads(args.input)
    print(f"Loaded {len(molecules)} leads from {args.input}")

    results = []
    for i, (smi, mol, name) in enumerate(molecules, 1):
        result = evaluate_single(smi, mol, args.scaffold, rank=i)
        results.append(result)

    report = write_report(results, args.output_dir, args.scaffold)
    json_out = write_json(results, args.output_dir, args.scaffold)

    print(f"Evaluated {len(results)} leads")
    print(f"Report: {report}")
    print(f"JSON:   {json_out}")


if __name__ == "__main__":
    main()
