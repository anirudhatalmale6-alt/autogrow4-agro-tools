"""
Descriptor engine for the Plant Systemic Bioavailability pipeline.

Computes all physicochemical descriptors needed by Gates 1-4, each
tagged with provenance (experimental / predicted_RDKit / predicted_SMARTS /
predicted_approx).  Experimental values can be injected via overrides so
the circularity guard (Section 5.3-5.4 of the Data Brief) is enforceable
in code.
"""

import math
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

# ── pKa SMARTS tables ──────────────────────────────────────────────────

ACIDIC_GROUPS = {
    "carboxylic_acid":  ("[CX3](=O)[O;H1]",              4.0),
    "sulfonic_acid":    ("[SX4](=O)(=O)[OH]",            -1.0),
    "phosphonic_acid1": ("[PX4](=O)([OH])([OH])",           2.2),
    "phosphonic_acid2": ("[PX4](=O)([OH])([O-])",         5.5),
    "tetrazole":        ("c1nnn[nH]1",                    4.7),
    "sulfonamide":      ("[NX3;H1][SX4](=O)=O",         10.0),
    "phenol":           ("c[OH]",                        10.0),
    "thiol":            ("[SH]",                          8.5),
}

BASIC_GROUPS = {
    "primary_amine":   ("[NX3;H2;!$(N-C=O)][CX4]",               10.5),
    "secondary_amine": ("[NX3;H1;!$(N-C=O)]([CX4])[CX4]",       10.7),
    "tertiary_amine":  ("[NX3;H0;!$(N-C=O)]([CX4])([CX4])[CX4]", 9.8),
    "pyridine":        ("c1ccncc1",                                5.2),
    "imidazole":       ("c1cnc[nH]1",                              6.0),
    "guanidine":       ("[NX3;!$(N[N+](=O)[O-]);!$(N~[N+])][CX3](=[NX2;!$(N~[N+])])[NX3;!$(N[N+](=O)[O-]);!$(N~[N+])]", 12.5),
}

# ── Labile-group SMARTS (controlled vocab, Gate 5 reactive-site tags) ──

LABILE_SMARTS = {
    "ester":               "[#6][CX3](=O)[OX2][#6]",
    "glucose_ester":       "[OX2][CX3](=O)[C@@H]1[C@H]([OH])[C@@H]([OH])[C@H]([OH])[C@@H]([OH])O1",
    "free_COOH":           "[CX3](=O)[OX2H1]",
    "free_phenol":         "c[OH]",
    "aliphatic_OH":        "[CX4][OH]",
    "aromatic_ring":       "c1ccccc1",
    "purine_N":            "c1ncnc2[nH]cnc12",
    "aromatic_amine":      "c[NH2]",
    "thiol":               "[#6][SH]",
    "epoxide":             "C1OC1",
    "michael_acceptor":    "[#6]=[#6][CX3]=O",
    "nitroaromatic":       "c[N+](=O)[O-]",
    "catechol":            "c([OH])c[OH]",
    "p_OH_cinnamate":      "[OH]c1ccc(/C=C/C(=O)[OH])cc1",
    "butenolide":          "O=C1C=COC1",
    "tertiary_amine":      "[NX3;H0;!$(N-C=O)]([#6])([#6])[#6]",
    "N_oxide":             "[NX4+][O-]",
    "phosphonate":         "[#6][PX4](=O)([OH])[OH]",
    "indole":              "c1ccc2[nH]ccc2c1",
    "quaternary_N":        "[NX4+;H0]",
    "steroidal_multi_OH":  "[C@@H]([OH])[C@H]([OH])",
    "even_C_acid_linker":  "[CX3](=O)[OH][CH2][CH2][CH2][CH2]",
}

# Pre-compile
_LABILE_COMPILED = {}


def _get_labile_patterns():
    if not _LABILE_COMPILED:
        for tag, sma in LABILE_SMARTS.items():
            pat = Chem.MolFromSmarts(sma)
            if pat is not None:
                _LABILE_COMPILED[tag] = pat
    return _LABILE_COMPILED


# ── pKa estimation ─────────────────────────────────────────────────────

def estimate_pka_groups(mol):
    """Return list of dicts with group name, type (acidic/basic), and pKa.

    Skips basic group matches where the key nitrogen already carries a
    permanent formal charge (e.g. N-methylpyridinium [n+]) — those are
    quaternary / permanently charged, not protonatable.
    """
    results = []
    for name, (sma, pka) in ACIDIC_GROUPS.items():
        pat = Chem.MolFromSmarts(sma)
        if pat and mol.HasSubstructMatch(pat):
            for i, match in enumerate(mol.GetSubstructMatches(pat)):
                results.append({"group": name, "type": "acidic",
                                "pKa_est": pka, "atoms": list(match),
                                "instance": i + 1})
                # Diprotic phosphonic acid: P(=O)(OH)(OH) has two pKas
                if name == "phosphonic_acid1":
                    results.append({"group": "phosphonic_acid2",
                                    "type": "acidic", "pKa_est": 5.5,
                                    "atoms": list(match), "instance": i + 1})
    for name, (sma, pka) in BASIC_GROUPS.items():
        pat = Chem.MolFromSmarts(sma)
        if pat and mol.HasSubstructMatch(pat):
            for i, match in enumerate(mol.GetSubstructMatches(pat)):
                n_atoms = [idx for idx in match
                           if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7]
                if any(mol.GetAtomWithIdx(idx).GetFormalCharge() != 0
                       for idx in n_atoms):
                    continue
                results.append({"group": name, "type": "basic",
                                "pKa_est": pka, "atoms": list(match),
                                "instance": i + 1})
    return results


# ── Net charge estimation ──────────────────────────────────────────────

def estimate_net_charge(mol, pH, pka_groups=None):
    """Estimate integer net formal charge at a given pH.

    Two contributions:
    1. pH-dependent ionisation from pKa groups (acid deprotonation / base
       protonation via Henderson-Hasselbalch).
    2. Permanent formal charges from the molecular graph (quaternary N+,
       metal centres, etc.) — these are pH-independent and always counted.
    """
    if pka_groups is None:
        pka_groups = estimate_pka_groups(mol)

    # atoms already accounted for by a pKa group
    pka_atom_set = set()
    for g in pka_groups:
        pka_atom_set.update(g["atoms"])

    charge = 0

    # pH-dependent ionisation
    for g in pka_groups:
        if g["type"] == "acidic" and pH > g["pKa_est"]:
            charge -= 1
        elif g["type"] == "basic" and pH < g["pKa_est"]:
            charge += 1

    # permanent formal charges (not covered by any pKa group)
    for atom in mol.GetAtoms():
        fc = atom.GetFormalCharge()
        if fc != 0 and atom.GetIdx() not in pka_atom_set:
            charge += fc

    return int(round(charge))


# ── Labile group detection ─────────────────────────────────────────────

def detect_labile_groups(mol):
    """Return sorted list of labile-group tags found in the molecule."""
    pats = _get_labile_patterns()
    found = []
    for tag, pat in pats.items():
        if mol.HasSubstructMatch(pat):
            found.append(tag)
    return sorted(found)


# ── ESOL aqueous solubility ────────────────────────────────────────────

def _esol_log_s(mol, logp=None):
    """ESOL log S in mol/L (Delaney 2004)."""
    if logp is None:
        logp = Crippen.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    rb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_heavy = mol.GetNumHeavyAtoms()
    ap = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic()) / n_heavy if n_heavy else 0
    return 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rb - 0.74 * ap


def _sol_mg_per_l(mol, log_s, mw=None):
    if mw is None:
        mw = Descriptors.MolWt(mol)
    return mw * (10 ** log_s) * 1000.0


# ── Abraham descriptor approximations ─────────────────────────────────

def _approx_abraham(mol):
    """Rough Abraham solute descriptors from RDKit properties.

    These are coarse placeholders.  Proper values should come from
    UFZ-LSER or experimental sources and be injected via overrides.
    """
    mr = Crippen.MolMR(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    mw = Descriptors.MolWt(mol)

    # E: excess molar refraction ≈ scaled MR
    E = (mr - 25.0) / 10.0

    # S: dipolarity/polarizability ≈ scaled TPSA
    S = min(tpsa / 60.0, 3.0)

    # A: H-bond acidity ≈ donor-count based
    A = min(hbd * 0.35, 2.5)

    # B: H-bond basicity ≈ acceptor-count based
    B = min(hba * 0.12, 2.5)

    # V: McGowan volume ≈ MW-based (McGowan & Abraham, 1987)
    n_heavy = mol.GetNumHeavyAtoms()
    n_bonds = mol.GetNumBonds()
    V = (mw / 100.0) - 0.0628 * n_bonds if n_heavy else 0

    return {"E": round(E, 3), "S": round(S, 3), "A": round(A, 3),
            "B": round(B, 3), "V": round(V, 3)}


# ── logD calculation ───────────────────────────────────────────────────

def _calc_log_d(logkow, pH, pka_acid=None, pka_base=None):
    """logD at a given pH from Henderson-Hasselbalch.

    For acids:  logD = logKow - log10(1 + 10^(pH - pKa))
    For bases:  logD = logKow - log10(1 + 10^(pKa - pH))
    For zwitterions both corrections apply additively.
    """
    correction = 0.0
    if pka_acid is not None:
        correction += math.log10(1.0 + 10 ** (pH - pka_acid))
    if pka_base is not None:
        correction += math.log10(1.0 + 10 ** (pka_base - pH))
    return logkow - correction


# ── Provenance flag ────────────────────────────────────────────────────

def descriptor_provenance_flag(desc_dict):
    """Summarise the provenance across all descriptors.

    Returns 'all_experimental', 'mixed', or 'all_predicted'.
    Ignores keys whose value is None (not computable / not applicable).
    """
    provs = set()
    for v in desc_dict.values():
        if isinstance(v, dict) and v.get("value") is not None:
            provs.add(v.get("provenance", "predicted"))
    if not provs:
        return "all_predicted"
    if provs == {"experimental"}:
        return "all_experimental"
    if "experimental" in provs:
        return "mixed"
    return "all_predicted"


# ── Main descriptor engine ─────────────────────────────────────────────

def compute_descriptors(mol, experimental_overrides=None):
    """Compute all descriptors for a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        An RDKit molecule object (already sanitised).
    experimental_overrides : dict, optional
        Keys are descriptor names, values are (value, provenance_string)
        tuples.  These replace the computed value.

    Returns
    -------
    dict[str, dict]
        Mapping descriptor_name → {"value": ..., "provenance": ...}.
    """
    ov = experimental_overrides or {}

    def _v(key, value, prov):
        if key in ov:
            val, p = ov[key]
            return {"value": val, "provenance": p}
        return {"value": value, "provenance": prov}

    # Basic RDKit descriptors
    logkow = Crippen.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    rb = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # pKa
    pka_groups = estimate_pka_groups(mol)
    acid_pkas = [g["pKa_est"] for g in pka_groups if g["type"] == "acidic"]
    base_pkas = [g["pKa_est"] for g in pka_groups if g["type"] == "basic"]
    pka_acid = min(acid_pkas) if acid_pkas else None
    pka_base = max(base_pkas) if base_pkas else None

    # logD at apoplast pH
    logd_55 = _calc_log_d(logkow, 5.5, pka_acid, pka_base)

    # Solubility
    log_s = _esol_log_s(mol, logkow)
    sol_mg = _sol_mg_per_l(mol, log_s, mw)

    # Net charges
    charge_75 = estimate_net_charge(mol, 7.5, pka_groups)
    charge_55 = estimate_net_charge(mol, 5.5, pka_groups)

    # Abraham approximations
    abr = _approx_abraham(mol)

    # Labile groups
    labile = detect_labile_groups(mol)

    d = {
        # Gate 1
        "logKow_value":      _v("logKow_value",      round(logkow, 2),   "predicted_RDKit"),
        "logD_pH5p5":        _v("logD_pH5p5",        round(logd_55, 2),  "predicted_RDKit"),
        "mw":                _v("mw",                round(mw, 1),       "predicted_RDKit"),
        "water_sol_mgL":     _v("water_sol_mgL",     round(sol_mg, 1),   "predicted_RDKit"),
        "melting_point_C":   _v("melting_point_C",   None,               None),
        "vapor_pressure_mPa": _v("vapor_pressure_mPa", None,             None),
        "abraham_E":         _v("abraham_E",         abr["E"],           "predicted_approx"),
        "abraham_S":         _v("abraham_S",         abr["S"],           "predicted_approx"),
        "abraham_A":         _v("abraham_A",         abr["A"],           "predicted_approx"),
        "abraham_B":         _v("abraham_B",         abr["B"],           "predicted_approx"),
        "abraham_V":         _v("abraham_V",         abr["V"],           "predicted_approx"),
        # Gate 2/3
        "pka_acidic":        _v("pka_acidic",        pka_acid,           "predicted_SMARTS" if pka_acid is not None else None),
        "pka_basic":         _v("pka_basic",         pka_base,           "predicted_SMARTS" if pka_base is not None else None),
        "tpsa":              _v("tpsa",              round(tpsa, 1),     "predicted_RDKit"),
        "hbd":               _v("hbd",               hbd,                "predicted_RDKit"),
        "hba":               _v("hba",               hba,                "predicted_RDKit"),
        "rotatable_bonds":   _v("rotatable_bonds",   rb,                 "predicted_RDKit"),
        "net_charge_pH7p5":  _v("net_charge_pH7p5",  charge_75,          "predicted_RDKit"),
        "net_charge_pH5p5":  _v("net_charge_pH5p5",  charge_55,          "predicted_RDKit"),
        # Gate 5
        "known_labile_groups": {"value": labile, "provenance": "predicted_SMARTS"},
    }

    # pKa group details (not a scored descriptor, but needed downstream)
    d["_pka_groups"] = {"value": pka_groups, "provenance": "predicted_SMARTS"}

    return d


# ── CLI test ───────────────────────────────────────────────────────────

SEED_COMPOUNDS = [
    ("PSB-0001", "Glyphosate",       "OC(=O)CNCP(=O)(O)O"),
    ("PSB-0002", "2,4-D",            "OC(=O)COc1ccc(Cl)cc1Cl"),
    ("PSB-0003", "IAA",              "OC(=O)Cc1c[nH]c2ccccc12"),
    ("PSB-0004", "Imidacloprid",     "O=N/[N]=C1\\NCCN1Cc1ccc(Cl)nc1"),
    ("PSB-0005", "Tebuconazole",     "CC(C)(C)C(O)(Cn1cncn1)CCc1ccc(Cl)cc1"),
    ("PSB-0006", "Nicotine",         "CN1CCC[C@H]1c1cccnc1"),
    ("PSB-0007", "Paraquat",         "C[n+]1ccc(-c2cc[n+](C)cc2)cc1"),
    ("PSB-0008", "Chlorothalonil",   "N#Cc1c(Cl)c(C#N)c(Cl)c(Cl)c1Cl"),
]


def _fmt(val, width=8):
    if val is None:
        return "---".rjust(width)
    if isinstance(val, float):
        return f"{val:>{width}.2f}"
    if isinstance(val, int):
        return f"{val:>{width}d}"
    if isinstance(val, list):
        s = ";".join(val) if val else "---"
        return s[:width*3].ljust(width*3)
    return str(val).rjust(width)


if __name__ == "__main__":
    print("=" * 120)
    print("PLANT SCORING PIPELINE - DESCRIPTOR ENGINE TEST")
    print("=" * 120)

    # Table 1: Gate 1 descriptors
    hdr = f"{'ID':<10} {'Name':<16} {'logKow':>7} {'logD55':>7} {'MW':>7} {'Sol_mg':>9} {'Abr_E':>7} {'Abr_S':>7} {'Abr_A':>7} {'Abr_B':>7} {'Abr_V':>7}"
    print("\nTABLE 1: GATE 1 - CUTICLE DESCRIPTORS")
    print("-" * len(hdr))
    print(hdr)
    print("-" * len(hdr))

    all_desc = []
    for cid, name, smi in SEED_COMPOUNDS:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"{cid:<10} {name:<16} *** PARSE FAILED ***")
            all_desc.append(None)
            continue
        desc = compute_descriptors(mol)
        all_desc.append(desc)
        g = lambda k: desc[k]["value"]
        print(f"{cid:<10} {name:<16} {_fmt(g('logKow_value'),7)} {_fmt(g('logD_pH5p5'),7)} "
              f"{_fmt(g('mw'),7)} {_fmt(g('water_sol_mgL'),9)} "
              f"{_fmt(g('abraham_E'),7)} {_fmt(g('abraham_S'),7)} {_fmt(g('abraham_A'),7)} "
              f"{_fmt(g('abraham_B'),7)} {_fmt(g('abraham_V'),7)}")

    # Table 2: Gate 2/3 descriptors
    print()
    hdr2 = f"{'ID':<10} {'Name':<16} {'pKa_a':>7} {'pKa_b':>7} {'TPSA':>7} {'HBD':>5} {'HBA':>5} {'RB':>4} {'Q@7.5':>6} {'Q@5.5':>6}"
    print("TABLE 2: GATE 2/3 - LOADING / TRANSPORT DESCRIPTORS")
    print("-" * len(hdr2))
    print(hdr2)
    print("-" * len(hdr2))

    for (cid, name, _), desc in zip(SEED_COMPOUNDS, all_desc):
        if desc is None:
            print(f"{cid:<10} {name:<16} *** SKIPPED ***")
            continue
        g = lambda k: desc[k]["value"]
        print(f"{cid:<10} {name:<16} {_fmt(g('pka_acidic'),7)} {_fmt(g('pka_basic'),7)} "
              f"{_fmt(g('tpsa'),7)} {g('hbd'):>5d} {g('hba'):>5d} {g('rotatable_bonds'):>4d} "
              f"{g('net_charge_pH7p5'):>6d} {g('net_charge_pH5p5'):>6d}")

    # Table 3: Labile groups
    print()
    print("TABLE 3: GATE 5 - LABILE GROUP TAGS")
    print("-" * 90)
    print(f"{'ID':<10} {'Name':<16} {'Labile Groups'}")
    print("-" * 90)

    for (cid, name, _), desc in zip(SEED_COMPOUNDS, all_desc):
        if desc is None:
            continue
        tags = desc["known_labile_groups"]["value"]
        print(f"{cid:<10} {name:<16} {'; '.join(tags) if tags else '---'}")

    # Provenance summary
    print()
    print("PROVENANCE SUMMARY")
    print("-" * 60)
    for (cid, name, _), desc in zip(SEED_COMPOUNDS, all_desc):
        if desc is None:
            continue
        flag = descriptor_provenance_flag(desc)
        print(f"  {cid} {name:<16} → {flag}")

    print()
    print("Done. All 8 seed compounds processed.")
