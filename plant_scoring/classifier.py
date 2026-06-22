"""
Module 7 — Compound Class Classifier (18 classes + GAP).

Routes each molecule to the correct Gate 3 mechanism and sets the
carrier-recognition probability (alpha) for Gate 2.  Detection is
SMARTS-first for structural classes, then descriptor-gated for
property-based classes.  Priority order ensures specific structural
classes outrank generic property classes.

Reference: Master Checkpoint v2, Section 4 + Section 5.
"""

import math
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

from plant_scoring.descriptors import (
    compute_descriptors, estimate_pka_groups, estimate_net_charge,
    _esol_log_s, _sol_mg_per_l,
)

# ── SMARTS library ────────────────────────────────────────────────────

_SMARTS_CACHE = {}

def _pat(name, smarts):
    if name not in _SMARTS_CACHE:
        p = Chem.MolFromSmarts(smarts)
        if p is None:
            raise ValueError(f"Bad SMARTS for {name}: {smarts}")
        _SMARTS_CACHE[name] = p
    return _SMARTS_CACHE[name]


# Structural SMARTS keyed to classes (compiled on first use)
STRUCTURAL_SMARTS = {
    # Phosphonate C-PO3H2
    "phosphonate":          "[#6][PX4](=O)([OH])[OH]",
    "phosphonate_mono":     "[#6][PX4](=O)([OH])[O-]",

    # Alpha-amino acid backbone: H2N-CH(R)-COOH
    "alpha_aa":             "[NX3;H2,H1][CX4H1]([*])[CX3](=O)[OX2H1,OX1-]",

    # Adenine core (purine with 6-amino / N6-substituted)
    "adenine_core":         "c1nc(N)c2nc[nH]c2n1",
    "adenine_n6sub":        "c1nc(N[#6])c2nc[nH]c2n1",
    # N6-substituted adenine with deprotonated N9 too
    "adenine_n6sub_v2":     "c1nc(N([#6])[#1,#6])c2ncnc2n1",

    # Phenylurea  NC(=O)Nc1ccccc1
    "phenylurea":           "[NX3][CX3](=O)[NX3]c1ccccc1",

    # Butenolide (gamma-butenolide / D-ring of strigolactones)
    "butenolide":           "O=C1C=COC1",

    # Strigolactone: tricyclic ABC + butenolide D
    # Simplified: fused bicyclic lactone connected via enol ether to butenolide
    "strigolactone_abc":    "C1CCC2C(C1)OC(=O)C2",

    # Iridoid: cyclopentane fused to pyran (iridoid skeleton)
    "cyclopentanoid_pyran": "C1CC2CCOC(O)C2C1",
    # beta-D-glucoside
    "glucoside":            "[OX2][C@@H]1[C@H]([OH])[C@@H]([OH])[C@H]([OH])[C@@H](CO)O1",
    # Simpler glucoside for non-stereo SMILES
    "glucoside_simple":     "OC1OC(CO)C(O)C(O)C1O",

    # Steroid ABCD ring system: four fused rings (6-6-6-5)
    "steroid_abcd":         "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~[#6]~2~[#6]~1",

    # Flavonoid scaffold: 2-phenylchromen-4-one
    "flavone":              "O=c1cc(-c2ccccc2)oc2ccccc12",
    "flavanone":            "O=C1CC(-c2ccccc2)Oc2ccccc21",
    "isoflavone":           "O=c1coc2ccccc2c1-c1ccccc1",
    "chalcone":             "O=C(/C=C/c1ccccc1)c1ccccc1",

    # Phenolic OH on aromatic ring
    "phenolic_oh":          "c[OH]",

    # Carboxylate (free carboxylic acid)
    "carboxylate":          "[CX3](=O)[OX2H1]",

    # Tertiary amine (not amide, not N-oxide)
    "tert_amine":           "[NX3;H0;!$(N-C=O);!$(N-[#16]);!$([N+])]([#6])([#6])[#6]",

    # N-oxide
    "n_oxide":              "[NX4+][O-]",

    # Bipyridinium / permanent poly-cation
    "n_plus_aromatic":      "[n+]",

    # Sesquiterpene: 15 carbons, 3 isoprene-derived rings — detect large
    # alicyclic system + methyl branches rather than exact skeleton
    "sesquiterpene_like":   "[CH3][C]1([CH3])CCC2CCCC([CH3])C2C1",

    # Diterpene: 20 carbons, polycyclic — broad pattern
    "diterpene_like":       "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6](~[#6]~[#6]~[#6]~3)~[#6]~2~[#6]~1",
}


# ── Class definitions ─────────────────────────────────────────────────

CLASS_PROPERTIES = {
    "PHOSPHONATE_NUTRIENT_MIMIC": {
        "carrier_alpha": 0.80,
        "gate3_mechanism": "carrier_override",
    },
    "AMINO_ACID_NEUTRAL": {
        "carrier_alpha": 0.80,
        "gate3_mechanism": "carrier_override",
    },
    "AMINO_ACID_ACIDIC": {
        "carrier_alpha": 0.80,
        "gate3_mechanism": "carrier_override",
    },
    "AMINO_ACID_BASIC": {
        "carrier_alpha": 0.70,
        "gate3_mechanism": "carrier_override",
    },
    "PURINE_CYTOKININ": {
        "carrier_alpha": 0.85,
        "gate3_mechanism": "carrier_override",
    },
    "PHENYLUREA_CYTOKININ": {
        "carrier_alpha": 0.15,
        "gate3_mechanism": "passive_moderate",
    },
    "STRIGOLACTONE_CANONICAL": {
        "carrier_alpha": -0.30,
        "gate3_mechanism": "efflux_penalty",
    },
    "IRIDOID_GLYCOSIDE": {
        "carrier_alpha": 0.65,
        "gate3_mechanism": "carrier_override",
    },
    "STEROID_HORMONE": {
        "carrier_alpha": 0.0,
        "gate3_mechanism": "negative_control",
    },
    "CARRIER_GA": {
        "carrier_alpha": 0.50,
        "gate3_mechanism": "partial_kleier",
    },
    "WEAK_ACID_SESQUITERPENOID": {
        "carrier_alpha": 0.20,
        "gate3_mechanism": "full_kleier",
    },
    "ALKALOID_N_OXIDE": {
        "carrier_alpha": 0.70,
        "gate3_mechanism": "carrier_override",
    },
    "PHENOLIC_FLAVONOID": {
        "carrier_alpha": 0.05,
        "gate3_mechanism": "partial_kleier",
    },
    "OTHER_permanent_dication": {
        "carrier_alpha": 0.0,
        "gate3_mechanism": "negative_control",
    },
    "TERTIARY_AMINE_ALKALOID": {
        "carrier_alpha": 0.0,
        "gate3_mechanism": "inverse_kleier",
    },
    "HIGH_LIPOPHILE": {
        "carrier_alpha": 0.0,
        "gate3_mechanism": "negative_control",
    },
    "ZWITTERION": {
        "carrier_alpha": 0.05,
        "gate3_mechanism": "full_kleier",
    },
    "WEAK_ACID_ORGANIC": {
        "carrier_alpha": 0.0,
        "gate3_mechanism": "full_kleier",
    },
    "POLAR_NEUTRAL_SYSTEMIC": {
        "carrier_alpha": 0.25,
        "gate3_mechanism": "intermediate_perm",
    },
    "NEUTRAL_LIPOPHILE": {
        "carrier_alpha": 0.0,
        "gate3_mechanism": "full_kleier",
    },
    "GAP": {
        "carrier_alpha": 0.0,
        "gate3_mechanism": "full_kleier",
    },
}


# ── Helper extractors ─────────────────────────────────────────────────

def _has(mol, name):
    return mol.HasSubstructMatch(_pat(name, STRUCTURAL_SMARTS[name]))

def _count(mol, name):
    return len(mol.GetSubstructMatches(_pat(name, STRUCTURAL_SMARTS[name])))

def _permanent_formal_charge(mol):
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())

def _ring_count(mol):
    return rdMolDescriptors.CalcNumRings(mol)

def _num_carbons(mol):
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)


def _extract_props(mol, descriptors):
    """Pull key numeric values from descriptor dict or compute on the fly."""
    if descriptors is not None:
        def _g(k):
            v = descriptors.get(k)
            if isinstance(v, dict):
                return v.get("value")
            return v
        logkow = _g("logKow_value")
        mw = _g("mw")
        pka_acid = _g("pka_acidic")
        pka_base = _g("pka_basic")
        sol = _g("water_sol_mgL")
        q75 = _g("net_charge_pH7p5")
        q55 = _g("net_charge_pH5p5")
    else:
        logkow = Crippen.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        pka_groups = estimate_pka_groups(mol)
        acids = [g["pKa_est"] for g in pka_groups if g["type"] == "acidic"]
        bases = [g["pKa_est"] for g in pka_groups if g["type"] == "basic"]
        pka_acid = min(acids) if acids else None
        pka_base = max(bases) if bases else None
        log_s = _esol_log_s(mol, logkow)
        sol = _sol_mg_per_l(mol, log_s, mw)
        q75 = estimate_net_charge(mol, 7.5, pka_groups)
        q55 = estimate_net_charge(mol, 5.5, pka_groups)

    return {
        "logkow": logkow if logkow is not None else Crippen.MolLogP(mol),
        "mw": mw if mw is not None else Descriptors.MolWt(mol),
        "pka_acid": pka_acid,
        "pka_base": pka_base,
        "sol": sol,
        "q75": q75 if q75 is not None else 0,
        "q55": q55 if q55 is not None else 0,
    }


# ── Classification cascade ────────────────────────────────────────────

def classify_compound(mol, descriptors=None):
    """Classify a molecule into one of the 18 Module 7 classes.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
    descriptors : dict, optional
        Output of compute_descriptors().  Computed internally if not given.

    Returns
    -------
    dict with module7_class, carrier_alpha, gate3_mechanism, confidence,
    matched_features.
    """
    p = _extract_props(mol, descriptors)

    # ── 1. PHOSPHONATE_NUTRIENT_MIMIC ──────────────────────────────
    if _has(mol, "phosphonate") or _has(mol, "phosphonate_mono"):
        if p["logkow"] < -1.0 or (p["pka_acid"] is not None and p["pka_acid"] < 6):
            return _result("PHOSPHONATE_NUTRIENT_MIMIC", 0.90,
                           ["C-PO3H2 group", f"logKow={p['logkow']:.1f}"])
        # Phosphonate but not very polar — still flag, lower confidence
        return _result("PHOSPHONATE_NUTRIENT_MIMIC", 0.60,
                       ["C-PO3H2 group", f"logKow={p['logkow']:.1f} (atypically high)"])

    # ── 2. AMINO_ACID_* ───────────────────────────────────────────
    # All amino acids have a backbone amine (pKa ~9-10) and backbone COOH
    # (pKa ~2).  The sub-class depends on the SIDE CHAIN character:
    # basic = extra basic group beyond backbone amine (Arg/Lys/His)
    # acidic = extra COOH beyond backbone (Glu/Asp)
    # neutral = everything else
    if _has(mol, "alpha_aa"):
        if descriptors is not None:
            pka_groups = descriptors.get("_pka_groups", {}).get("value", [])
        else:
            pka_groups = estimate_pka_groups(mol)
        n_acid = sum(1 for g in pka_groups if g["type"] == "acidic")
        n_base = sum(1 for g in pka_groups if g["type"] == "basic")
        # Extra basic group beyond backbone primary amine
        has_sidechain_base = n_base >= 2 or (
            n_base >= 1 and any(g["group"] in ("guanidine", "imidazole")
                                for g in pka_groups if g["type"] == "basic"))
        # Extra acidic group beyond backbone COOH
        has_sidechain_acid = n_acid >= 2
        if has_sidechain_base:
            return _result("AMINO_ACID_BASIC", 0.85,
                           ["alpha-AA backbone", "basic side chain",
                            f"{n_base} basic groups"])
        if has_sidechain_acid:
            return _result("AMINO_ACID_ACIDIC", 0.85,
                           ["alpha-AA backbone", "acidic side chain",
                            f"{n_acid} acidic groups"])
        return _result("AMINO_ACID_NEUTRAL", 0.80,
                       ["alpha-AA backbone"])

    # ── 3. PURINE_CYTOKININ ───────────────────────────────────────
    if _has(mol, "adenine_n6sub") or _has(mol, "adenine_n6sub_v2"):
        return _result("PURINE_CYTOKININ", 0.90,
                       ["adenine core", "N6-substitution"])
    if _has(mol, "adenine_core") and p["mw"] > 200:
        return _result("PURINE_CYTOKININ", 0.65,
                       ["adenine core", f"MW={p['mw']:.0f}"])

    # ── 4. PHENYLUREA_CYTOKININ ───────────────────────────────────
    if _has(mol, "phenylurea"):
        return _result("PHENYLUREA_CYTOKININ", 0.80,
                       ["phenylurea moiety"])

    # ── 5. STRIGOLACTONE_CANONICAL ────────────────────────────────
    if _has(mol, "butenolide"):
        if _has(mol, "strigolactone_abc") or _ring_count(mol) >= 3:
            return _result("STRIGOLACTONE_CANONICAL", 0.80,
                           ["butenolide D-ring", f"{_ring_count(mol)} rings"])
        # Butenolide alone without the ABC system — partial match
        if _ring_count(mol) >= 2 and _num_carbons(mol) >= 15:
            return _result("STRIGOLACTONE_CANONICAL", 0.50,
                           ["butenolide D-ring", "polycyclic terpene-like"])

    # ── 6. IRIDOID_GLYCOSIDE ─────────────────────────────────────
    has_glucoside = _has(mol, "glucoside") or _has(mol, "glucoside_simple")
    if has_glucoside and _has(mol, "cyclopentanoid_pyran"):
        return _result("IRIDOID_GLYCOSIDE", 0.85,
                       ["cyclopentanoid-pyran", "beta-D-glucoside"])
    # Looser: glucoside + cyclopentane ring + MW in iridoid range
    if has_glucoside and _ring_count(mol) >= 2 and 300 < p["mw"] < 500:
        n_5ring = sum(1 for r in mol.GetRingInfo().AtomRings() if len(r) == 5)
        if n_5ring >= 1:
            return _result("IRIDOID_GLYCOSIDE", 0.50,
                           ["beta-D-glucoside", "5-membered ring", f"MW={p['mw']:.0f}"])

    # ── 7. STEROID_HORMONE ────────────────────────────────────────
    if _has(mol, "steroid_abcd"):
        if p["pka_acid"] is None and p["pka_base"] is None:
            return _result("STEROID_HORMONE", 0.85,
                           ["steroidal ABCD ring", "no ionizable group"])
        # Steroid with ionizable group — could be a steroidal acid
        return _result("STEROID_HORMONE", 0.55,
                       ["steroidal ABCD ring", "has ionizable group (atypical)"])

    # ── 8. CARRIER_GA (diterpene + COOH, MW 300-450) ─────────────
    if (_has(mol, "diterpene_like") and _has(mol, "carboxylate")
            and 300 <= p["mw"] <= 450):
        return _result("CARRIER_GA", 0.75,
                       ["diterpene skeleton", "free COOH", f"MW={p['mw']:.0f}"])

    # ── 9. WEAK_ACID_SESQUITERPENOID ─────────────────────────────
    if _has(mol, "sesquiterpene_like") and _has(mol, "carboxylate"):
        return _result("WEAK_ACID_SESQUITERPENOID", 0.70,
                       ["sesquiterpene skeleton", "free COOH"])
    # ABA-like: ~15 carbons, cyclohexenone, COOH
    if (_has(mol, "carboxylate") and 12 <= _num_carbons(mol) <= 18
            and _ring_count(mol) >= 1 and p["pka_acid"] is not None
            and 3.5 <= p["pka_acid"] <= 6.0):
        n_5_6 = sum(1 for r in mol.GetRingInfo().AtomRings() if len(r) in (5, 6))
        if n_5_6 >= 1 and p["mw"] < 350:
            return _result("WEAK_ACID_SESQUITERPENOID", 0.40,
                           ["COOH + terpene-like C count", f"C={_num_carbons(mol)}"])

    # ── 10. ALKALOID_N_OXIDE ─────────────────────────────────────
    if _has(mol, "n_oxide"):
        return _result("ALKALOID_N_OXIDE", 0.85,
                       ["R3N+-O- group"])

    # ── 11. PHENOLIC_FLAVONOID ───────────────────────────────────
    is_flavonoid = any(_has(mol, k) for k in
                       ("flavone", "flavanone", "isoflavone", "chalcone"))
    if is_flavonoid and _has(mol, "phenolic_oh"):
        return _result("PHENOLIC_FLAVONOID", 0.85,
                       ["flavonoid scaffold", "phenolic OH"])
    if is_flavonoid:
        return _result("PHENOLIC_FLAVONOID", 0.70,
                       ["flavonoid scaffold"])

    # ── 12. OTHER_permanent_dication ─────────────────────────────
    pfc = _permanent_formal_charge(mol)
    if pfc >= 2:
        n_plus = _count(mol, "n_plus_aromatic")
        feats = [f"permanent charge={pfc:+d}"]
        if n_plus >= 2:
            feats.append(f"{n_plus}x aromatic N+")
        return _result("OTHER_permanent_dication", 0.90, feats)

    # ── 13. TERTIARY_AMINE_ALKALOID ──────────────────────────────
    if _has(mol, "tert_amine"):
        if p["pka_base"] is not None and 6.0 <= p["pka_base"] <= 10.0:
            return _result("TERTIARY_AMINE_ALKALOID", 0.80,
                           ["tertiary amine", f"pKa_base={p['pka_base']:.1f}"])
    # Pyridine/imidazole base in alkaloid range
    if p["pka_base"] is not None and 6.0 <= p["pka_base"] <= 10.0:
        if p["pka_acid"] is None and _ring_count(mol) >= 2:
            return _result("TERTIARY_AMINE_ALKALOID", 0.55,
                           ["heterocyclic base", f"pKa_base={p['pka_base']:.1f}"])

    # ── 14. HIGH_LIPOPHILE ───────────────────────────────────────
    if p["logkow"] > 5.0:
        return _result("HIGH_LIPOPHILE", 0.85,
                       [f"logKow={p['logkow']:.1f} (>5)"])

    # ── 15. ZWITTERION ───────────────────────────────────────────
    if (p["pka_acid"] is not None and 2.0 <= p["pka_acid"] <= 7.0
            and p["pka_base"] is not None and 5.0 <= p["pka_base"] <= 10.0):
        return _result("ZWITTERION", 0.70,
                       [f"acid pKa={p['pka_acid']:.1f}", f"base pKa={p['pka_base']:.1f}"])

    # ── 16. WEAK_ACID_ORGANIC ────────────────────────────────────
    if p["pka_acid"] is not None and 2.0 <= p["pka_acid"] <= 7.0:
        feats = [f"acidic pKa={p['pka_acid']:.1f}"]
        if _has(mol, "carboxylate"):
            feats.append("free COOH")
        elif _has(mol, "phenolic_oh"):
            feats.append("phenol")
        return _result("WEAK_ACID_ORGANIC", 0.80, feats)

    # ── 17. POLAR_NEUTRAL_SYSTEMIC ───────────────────────────────
    # Effectively neutral at sieve-tube pH: no acid, and either no base
    # or only a weak base (pKa < 6) that is mostly unprotonated at pH 7.5.
    effectively_neutral = (
        p["pka_acid"] is None
        and (p["pka_base"] is None or p["pka_base"] < 6.0)
    )
    # logKow window: spec says -0.5 to 0.5, but imidacloprid (the defining
    # example, logKow 0.57) sits right at the edge.  Use -0.5 to 1.0 and
    # let the calibration step tighten if needed.
    if effectively_neutral and -0.5 <= p["logkow"] <= 1.0:
        sol = p["sol"]
        if sol is not None and sol > 500:
            feats = [f"logKow={p['logkow']:.2f}", f"sol={sol:.0f} mg/L"]
            if p["pka_base"] is not None:
                feats.append(f"weak base pKa={p['pka_base']:.1f} (effectively neutral)")
            else:
                feats.append("neutral")
            return _result("POLAR_NEUTRAL_SYSTEMIC", 0.75, feats)
    if effectively_neutral and -0.5 <= p["logkow"] <= 1.0 and p["mw"] < 350:
        sol = p["sol"]
        if sol is not None and sol > 500:
            return _result("POLAR_NEUTRAL_SYSTEMIC", 0.55,
                           [f"logKow={p['logkow']:.2f}", f"sol={sol:.0f} mg/L",
                            "near-neutral"])

    # ── 18. NEUTRAL_LIPOPHILE ────────────────────────────────────
    if p["logkow"] > 2.5:
        if p["pka_acid"] is None:
            return _result("NEUTRAL_LIPOPHILE", 0.75,
                           [f"logKow={p['logkow']:.1f}", "no ionizable acid"])
        # Has weak ionizable group but high logKow dominates
        if p["pka_acid"] > 7.0:
            return _result("NEUTRAL_LIPOPHILE", 0.55,
                           [f"logKow={p['logkow']:.1f}", f"weak acid pKa={p['pka_acid']:.1f} (>7)"])

    # ── 19. GAP (fallback) ───────────────────────────────────────
    feats = []
    if p["logkow"] is not None:
        feats.append(f"logKow={p['logkow']:.1f}")
    if p["pka_acid"] is not None:
        feats.append(f"pKa_acid={p['pka_acid']:.1f}")
    if p["pka_base"] is not None:
        feats.append(f"pKa_base={p['pka_base']:.1f}")
    feats.append("no class matched")
    return _result("GAP", 0.10, feats)


def _result(cls, confidence, features):
    props = CLASS_PROPERTIES[cls]
    return {
        "module7_class": cls,
        "carrier_alpha": props["carrier_alpha"],
        "gate3_mechanism": props["gate3_mechanism"],
        "confidence": round(confidence, 2),
        "matched_features": features,
    }


def reclassify_after_transform(product_mol, descriptors=None):
    """Re-classify a product molecule after a Gate 5 class-switching transform.

    Used after B11 (sulfation → permanent anion) or E1 (nitro-reduction
    → amine), which change the compound's Gate 3 mechanism class
    mid-cascade.
    """
    return classify_compound(product_mol, descriptors)


# ── CLI test ──────────────────────────────────────────────────────────

SEED_COMPOUNDS = [
    ("PSB-0001", "Glyphosate",       "OC(=O)CNCP(=O)(O)O",
     "PHOSPHONATE_NUTRIENT_MIMIC"),
    ("PSB-0002", "2,4-D",            "OC(=O)COc1ccc(Cl)cc1Cl",
     "WEAK_ACID_ORGANIC"),
    ("PSB-0003", "IAA",              "OC(=O)Cc1c[nH]c2ccccc12",
     "WEAK_ACID_ORGANIC"),
    ("PSB-0004", "Imidacloprid",     "O=[N+]([O-])NC1=NCCN1Cc1ccc(Cl)nc1",
     "POLAR_NEUTRAL_SYSTEMIC"),
    ("PSB-0005", "Tebuconazole",     "CC(C)(C)C(O)(Cn1cncn1)CCc1ccc(Cl)cc1",
     "NEUTRAL_LIPOPHILE"),
    ("PSB-0006", "Nicotine",         "CN1CCC[C@H]1c1cccnc1",
     "TERTIARY_AMINE_ALKALOID"),
    ("PSB-0007", "Paraquat",         "C[n+]1ccc(-c2cc[n+](C)cc2)cc1",
     "OTHER_permanent_dication"),
    ("PSB-0008", "Chlorothalonil",   "N#Cc1c(Cl)c(C#N)c(Cl)c(Cl)c1Cl",
     "NEUTRAL_LIPOPHILE"),
]


if __name__ == "__main__":
    print("=" * 110)
    print("MODULE 7 — COMPOUND CLASS CLASSIFIER TEST")
    print("=" * 110)
    print(f"\n{'ID':<10} {'Name':<16} {'Expected':<30} {'Predicted':<30} {'Match':>5}  {'Conf':>4}  Gate3 Mechanism")
    print("-" * 110)

    correct = 0
    total = len(SEED_COMPOUNDS)

    for cid, name, smi, expected in SEED_COMPOUNDS:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"{cid:<10} {name:<16} *** PARSE FAILED ***")
            continue

        desc = compute_descriptors(mol)
        result = classify_compound(mol, desc)
        predicted = result["module7_class"]
        match = predicted == expected
        if match:
            correct += 1
        marker = "OK" if match else "MISS"

        print(f"{cid:<10} {name:<16} {expected:<30} {predicted:<30} {marker:>5}  "
              f"{result['confidence']:.2f}  {result['gate3_mechanism']}")
        if result["matched_features"]:
            print(f"{'':>10} {'':>16}   features: {', '.join(result['matched_features'])}")

    print(f"\nAccuracy: {correct}/{total} ({100*correct/total:.0f}%)")

    # Also test a class-switching scenario: sulfated salicylate
    print("\n--- Class-switching test: sulfated salicylic acid ---")
    # SA sulfate ester: installs permanent anion
    sa_sulfate = Chem.MolFromSmiles("O=C(O)c1ccccc1OS(=O)(=O)O")
    if sa_sulfate:
        r = classify_compound(sa_sulfate)
        print(f"  Salicylate-sulfate → {r['module7_class']} "
              f"(alpha={r['carrier_alpha']:.2f}, {r['gate3_mechanism']})")
        print(f"  features: {', '.join(r['matched_features'])}")
