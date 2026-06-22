"""
Gate 5 — Metabolic Stability (Module 6 v2)

Reactive-site alert library with SMARTS detection and the penalty multiplier
model from the Master Checkpoint v2, Section 6.

penalty(alert) = severity × probability × (1 − recovery)
Gate5_multiplier = ∏(1 − penalty_i), bounded 0–1
ACTIVATING alerts contribute zero penalty.
"""

from rdkit import Chem

HIGH = 0.8
MED = 0.5
LOW = 0.2

# ── Alert Library ─────────────────────────────────────────────────────
#
# Each entry: smarts, severity, probability, recovery, direction,
#             gate3_class_switch, description
#
# Direction key:
#   DEG = degradation/inactivation (penalised)
#   INACT = inactivation to vacuole (penalised, often reversible)
#   ACT / ACTIVATE = prodrug activation (zero penalty)
#   BUFFER = reversible storage (zero or negligible penalty)
#   STABLE = compound resists this route (zero penalty)
#   FORM_SWITCH = reversible interconversion (very low penalty)
#   TRANSFORM = product feeds next phase (penalised at reduced rate)

ALERT_LIBRARY = {
    # ── Class A: Oxidative (Phase I) ──────────────────────────────────

    "A1_indole": {
        "smarts": "c1ccc2[nH]ccc2c1",
        "severity": HIGH, "probability": HIGH, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Indole oxidation → OxIAA / phenoxy radical (peroxidase)",
    },
    "A1_phenol": {
        "smarts": "[OH]c1cc[cH]cc1",
        "severity": HIGH, "probability": HIGH, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Oxidisable phenol → phenoxy radical (class III Prx)",
    },
    "A2": {
        "smarts": "[nH]1c2ncnc2nc1[CH2,CH]C=C",
        "severity": HIGH, "probability": HIGH, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Adenine-N6 + isoprenyl C=C → CKX side-chain cleavage",
    },
    "A3": {
        "smarts": "O=C1C=CC(C)(C(/O)=C/C(=O)O)CC1",
        "severity": HIGH, "probability": HIGH, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Sesquiterpene cyclohexenone (ABA-like) → 8'-OH (CYP707A)",
    },
    "A4": {
        "smarts": "[C@@H]1CC[C@@H]2[C@H]1CC[C@H]1[C@H]2CC[C@@H]2C[C@@H](O)CC[C@]12C",
        "severity": HIGH, "probability": HIGH, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Steroidal ABCD ring → C-26 OH (CYP734A1/BAS1)",
    },
    "A5": {
        "smarts": "O=C1CCC[C@@H]1C/C=C",
        "severity": HIGH, "probability": MED, "recovery": 0.0,
        "direction": "INACT",
        "gate3_class_switch": False,
        "description": "Cyclopentanone oxylipin → ω-C12 hydroxylation",
    },
    "A7": {
        "smarts": "c([OH])c[OH]",
        "severity": MED, "probability": MED, "recovery": 0.0,
        "direction": "PROD_MOBILE",
        "gate3_class_switch": False,
        "description": "Catechol → O-methylation (COMT) → ferulic (mobile product)",
    },
    "A8": {
        "smarts": "[CH3]C(=C)CCC=C(C)C",
        "severity": MED, "probability": MED, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Volatile monoterpene → CYP71 hydroxylation + volatility loss",
    },
    "A9": {
        "smarts": "[#7;R;!$([#7+])]([CH3,CH2])~[#6;R]",
        "severity": MED, "probability": MED, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Heterocyclic N-alkyl → ring-OH / N-dealkylation (CYP450)",
    },
    "A10": {
        "smarts": "OC(=C(O)C(=O)CO)C(O)=O",
        "severity": MED, "probability": HIGH, "recovery": 0.5,
        "direction": "BUFFER",
        "gate3_class_switch": False,
        "description": "Ene-diol (ascorbate) → MDHA/DHA (reversible redox buffer)",
    },
    "A11": {
        "smarts": "[OH]c1ccc(/C=C/C(=O)O)cc1",
        "severity": HIGH, "probability": HIGH, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "p-OH-cinnamate/monolignol → cell-wall bound residue (Prx)",
    },

    # ── Class B: Conjugative (Phase II) ───────────────────────────────

    "B1": {
        "smarts": "c[OH]",
        "severity": MED, "probability": HIGH, "recovery": 0.5,
        "direction": "INACT",
        "gate3_class_switch": False,
        "description": "Aromatic -OH → O-glucoside (UGT D/E, ~44 enzymes)",
    },
    "B2": {
        "smarts": "[CX4][OH]",
        "severity": MED, "probability": MED, "recovery": 0.5,
        "direction": "INACT",
        "gate3_class_switch": False,
        "description": "Aliphatic -OH → O-glucoside (UGT)",
    },
    "B3": {
        "smarts": "[CX3](=O)[OX2H1]",
        "severity": MED, "probability": MED, "recovery": 0.5,
        "direction": "INACT",
        "gate3_class_switch": False,
        "description": "Free -COOH → glucose ester (UGT74F2/71B-C)",
    },
    "B4": {
        "smarts": "c1ncnc2[nH]cnc12",
        "severity": HIGH, "probability": MED, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Purine N7/N9 → N-glucoside (UGT76C1/C2, irreversible)",
    },
    "B5": {
        "smarts": "c[NH2]",
        "severity": HIGH, "probability": LOW, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Aromatic amine (aniline) → N-glucoside (UGT72B1 only)",
    },
    "B6": {
        "smarts": "[#6][SH]",
        "severity": MED, "probability": MED, "recovery": 0.25,
        "direction": "INACT",
        "gate3_class_switch": False,
        "description": "Free -SH (thiol) → S-glucoside (UGT D/E)",
    },
    "B7a": {
        "smarts": "[C;R]1[O;R][C;R]1",
        "severity": HIGH, "probability": HIGH, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Epoxide / arene oxide → GSH thioether (GST τ/φ, irreversible)",
    },
    "B7a_halide": {
        "smarts": "[CX4][F,Cl,Br,I]",
        "severity": HIGH, "probability": HIGH, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Activated halide → GSH displacement (GST τ/φ)",
    },
    "B7b": {
        "smarts": "[#6]=[#6][CX3]=O",
        "severity": MED, "probability": HIGH, "recovery": 0.5,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "α,β-unsat carbonyl / Michael acceptor → GSH β-thioether (retro-Michael reversible)",
    },
    "B8": {
        "smarts": "[C@@H]([OH])[C@H]([OH])",
        "severity": MED, "probability": MED, "recovery": 0.5,
        "direction": "INACT",
        "gate3_class_switch": False,
        "description": "Steroidal multi-OH → positional O-glucosides",
    },
    "B9": {
        "smarts": "[OH]c1ccc(C=CC(=O)[OH])cc1",
        "severity": HIGH, "probability": MED, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "p-OH cinnamic + free -COOH → CoA thioester → lignin (4CL)",
    },
    "B10_IAA": {
        "smarts": "OC(=O)Cc1c[nH]c2ccccc12",
        "severity": HIGH, "probability": HIGH, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Indole-acetic acid -COOH → GH3 amide conjugate (DEG for IAA)",
    },
    "B10_JA": {
        "smarts": "O=C1CCC[C@@H]1CC",
        "severity": 0.0, "probability": 0.0, "recovery": 0.0,
        "direction": "ACTIVATE",
        "gate3_class_switch": False,
        "description": "Cyclopentanone oxylipin -COOH → GH3 amide (ACT for JA)",
    },
    "B11": {
        "smarts": "c[OH]",
        "severity": 0.65, "probability": MED, "recovery": 0.0,
        "direction": "INACT",
        "gate3_class_switch": True,
        "description": "Phenolic OH → sulfate ester (SOT, PAPS) → class switch to permanent anion",
    },
    "B12": {
        "smarts": "c[OH]",
        "severity": MED, "probability": MED, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "β-D-glucoside 6'-OH → malonyl-glucoside (BAHD) — vacuole target; fires with B1/B2",
    },

    # ── Class C: Hydrolytic ───────────────────────────────────────────

    "C1": {
        "smarts": "c1ccc2[nH]ccc2c1CCCC",
        "severity": 0.0, "probability": 0.0, "recovery": 0.0,
        "direction": "ACTIVATE",
        "gate3_class_switch": False,
        "description": "Indole + 4-C side chain → β-oxidation → IAA (prodrug activation)",
    },
    "C2": {
        "smarts": "O=C(O/C=C/c1cc(O)c(O)cc1)C1CC(O)(C(=O)O)CC(O)C1O",
        "severity": MED, "probability": MED, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Caffeoyl-quinic ester (chlorogenic acid) → ester hydrolysis",
    },
    "C6": {
        "smarts": "[#6][CX3](=O)[OX2][CH2,CH;!$([OH])]",
        "severity": MED, "probability": MED, "recovery": 0.0,
        "direction": "CONTEXT",
        "gate3_class_switch": False,
        "description": "Ester R-C(=O)-O-R' → acid + alcohol (CXE carboxylesterase)",
    },
    "C7": {
        "smarts": "[CX3](=O)[OH]CCCC",
        "severity": 0.0, "probability": 0.0, "recovery": 0.0,
        "direction": "ACTIVATE",
        "gate3_class_switch": False,
        "description": "Terminal -COOH + even-C linker → β-oxidation chain-shortening (prodrug)",
    },

    # ── Class D: Non-enzymatic / form-switch ──────────────────────────

    "D1": {
        "smarts": "O=C1C=COC1",
        "severity": HIGH, "probability": HIGH, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Butenolide D-ring → pH-dependent Michael add.-elim. (maximal at phloem pH)",
    },
    "D2": {
        "smarts": "[CX3](=O)[OX2][CH3]",
        "severity": LOW, "probability": MED, "recovery": 0.5,
        "direction": "FORM_SWITCH",
        "gate3_class_switch": False,
        "description": "Ar/oxylipin methyl ester ⇄ free acid (SABATH fwd / MES rev)",
    },
    "D3": {
        "smarts": "OC(=O)c1ccccc1O",
        "severity": HIGH, "probability": MED, "recovery": 0.0,
        "direction": "DEG",
        "gate3_class_switch": False,
        "description": "Salicylate → 3-hydroxylation → 2,3-DHBA (S3H, senescence)",
    },
    "D4": {
        "smarts": "[NX3;H0;!$(N-C=O)]([#6])([#6])[#6]",
        "severity": LOW, "probability": LOW, "recovery": 0.0,
        "direction": "STABLE",
        "gate3_class_switch": False,
        "description": "Tertiary-amine alkaloid — minimal plant catabolism (stable)",
    },
    "D5": {
        "smarts": "[NX4+][O-]",
        "severity": 0.0, "probability": 0.0, "recovery": 0.0,
        "direction": "ACTIVATE",
        "gate3_class_switch": False,
        "description": "N-oxide → N-oxide reductase in sink (re-activation, not catabolism)",
    },
    "D6": {
        "smarts": "[#6][PX4](=O)([OH])[OH]",
        "severity": LOW, "probability": LOW, "recovery": 0.0,
        "direction": "STABLE",
        "gate3_class_switch": False,
        "description": "Phosphonate C-PO₃H₂ — no C-P lyase in higher plants (stable)",
    },

    # ── Class E: Reductive ────────────────────────────────────────────

    "E1": {
        "smarts": "c[N+](=O)[O-]",
        "severity": MED, "probability": 0.35,
        "recovery": 0.0,
        "direction": "TRANSFORM",
        "gate3_class_switch": True,
        "description": "Ar-NO₂ → Ar-NH₂ (nitroreductase/OYE) → Phase II handoff; class switch",
    },
    "E2": {
        "smarts": "[#6]=[#6]C(=O)[#6]",
        "severity": LOW, "probability": MED, "recovery": 0.0,
        "direction": "DETOXIFY",
        "gate3_class_switch": False,
        "description": "α,β-unsaturated carbonyl → saturated (OYE/OPR); competes with B7b",
    },
}

# ── Pre-compiled patterns ─────────────────────────────────────────────

_COMPILED = {}


def _compile():
    if _COMPILED:
        return _COMPILED
    for aid, spec in ALERT_LIBRARY.items():
        pat = Chem.MolFromSmarts(spec["smarts"])
        if pat is not None:
            _COMPILED[aid] = pat
        else:
            pass  # invalid SMARTS silently skipped; logged below in detect
    return _COMPILED


# ── Zero-penalty directions ───────────────────────────────────────────

_ZERO_PENALTY_DIRS = {"ACTIVATE", "BUFFER", "STABLE", "FORM_SWITCH"}


def _is_activating(spec):
    if spec["direction"] in _ZERO_PENALTY_DIRS:
        return True
    if spec["severity"] == 0.0:
        return True
    return False


# ── Public API ────────────────────────────────────────────────────────

def get_alert_info(alert_id):
    return ALERT_LIBRARY.get(alert_id)


def detect_alerts(mol):
    """Detect all firing Gate 5 alerts in a molecule.

    Returns list of dicts, each with:
        alert_id, smarts, matched_atoms (first match), direction,
        severity, probability, recovery, gate3_class_switch, description
    """
    pats = _compile()
    fired = []
    for aid, pat in pats.items():
        spec = ALERT_LIBRARY[aid]
        if mol.HasSubstructMatch(pat):
            match = mol.GetSubstructMatch(pat)
            fired.append({
                "alert_id": aid,
                "smarts": spec["smarts"],
                "matched_atoms": list(match),
                "direction": spec["direction"],
                "severity": spec["severity"],
                "probability": spec["probability"],
                "recovery": spec["recovery"],
                "gate3_class_switch": spec["gate3_class_switch"],
                "description": spec["description"],
            })
    return fired


def compute_gate5_multiplier(mol, alerts=None, classification=None):
    """Compute the Gate 5 metabolic stability multiplier.

    Parameters
    ----------
    mol : RDKit Mol
    alerts : list of dicts from detect_alerts(), or None (computed)
    classification : dict from classify_compound() (optional, used for
        B10 branch logic and B12 co-fire with B1/B2)

    Returns
    -------
    (multiplier, details) where details is a dict with fired_alerts,
    penalties, activating_alerts, class_switch_alerts, total_penalty_load
    """
    if alerts is None:
        alerts = detect_alerts(mol)

    fired_ids = [a["alert_id"] for a in alerts]

    # B12 only fires when B1 or B2 also fires
    b1_fires = "B1" in fired_ids
    b2_fires = "B2" in fired_ids
    b12_fires = "B12" in fired_ids

    penalties = {}
    activating = []
    class_switch = []

    for a in alerts:
        aid = a["alert_id"]
        spec = ALERT_LIBRARY[aid]

        if _is_activating(spec):
            activating.append(aid)
            penalties[aid] = 0.0
            continue

        sev = spec["severity"]
        prob = spec["probability"]
        rec = spec["recovery"]

        # B12 only counts when B1 or B2 co-fire
        if aid == "B12":
            if not (b1_fires or b2_fires):
                activating.append(aid)
                penalties[aid] = 0.0
                continue

        # When B12 co-fires, reduce B1/B2 recovery from 0.5 to 0.15
        if aid in ("B1", "B2") and b12_fires:
            rec = 0.15

        # A7 catechol → ferulic is beneficial; reduce penalty heavily
        if aid == "A7":
            sev = LOW
            prob = LOW

        # A1_phenol: don't double-count if A11 already fires on the
        # same p-OH-cinnamate scaffold
        if aid == "A1_phenol" and "A11" in fired_ids:
            penalties[aid] = 0.0
            continue

        # B11 sulfation: don't double-count with B1 on the same phenol
        if aid == "B11" and "B1" in fired_ids:
            # B11 fires independently but both can't fully penalise
            # the same pool; halve B11's probability
            prob = prob * 0.5

        # E2 competes with B7b for the same Michael acceptor warhead
        if aid == "E2" and "B7b" in fired_ids:
            prob = prob * 0.5

        pen = sev * prob * (1.0 - rec)
        penalties[aid] = round(pen, 4)

        if spec["gate3_class_switch"]:
            class_switch.append(aid)

    # Multiplier
    multiplier = 1.0
    for pen in penalties.values():
        multiplier *= (1.0 - pen)
    multiplier = max(0.0, min(1.0, multiplier))

    total_load = sum(penalties.values())

    return round(multiplier, 4), {
        "fired_alerts": fired_ids,
        "penalties": penalties,
        "activating_alerts": activating,
        "class_switch_alerts": class_switch,
        "total_penalty_load": round(total_load, 4),
    }


# ── Demo / test ───────────────────────────────────────────────────────

if __name__ == "__main__":
    TESTS = [
        ("OC(=O)CNCP(=O)(O)O", "Glyphosate"),
        ("OC(=O)COc1ccc(Cl)cc1Cl", "2,4-D"),
        ("OC(=O)Cc1c[nH]c2ccccc12", "IAA"),
        ("O=[N+]([O-])NC1=NCCN1Cc1ccc(Cl)nc1", "Imidacloprid"),
        ("CC(C)(C)C(O)(CCc1ccc(Cl)cc1)Cn1cncn1", "Tebuconazole"),
        ("CN1CCC[C@H]1c1cccnc1", "Nicotine"),
        ("C[n+]1ccc(-c2cc[n+](C)cc2)cc1", "Paraquat"),
        ("N#Cc1c(Cl)c(C#N)c(Cl)c(Cl)c1Cl", "Chlorothalonil"),
        ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin"),
        ("O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12", "Quercetin"),
        ("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "Caffeine"),
        ("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]", "TNT"),
    ]

    print(f"{'Name':<16} {'Mult':>5} {'#Fire':>5} {'#Act':>4} {'#ClSw':>5} "
          f"{'Penalty':>7}  Fired alerts")
    print("=" * 100)
    for smi, name in TESTS:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"{name:<16} INVALID SMILES")
            continue
        mult, det = compute_gate5_multiplier(mol)
        n_fire = len(det["fired_alerts"])
        n_act = len(det["activating_alerts"])
        n_cs = len(det["class_switch_alerts"])
        load = det["total_penalty_load"]
        aids = ", ".join(det["fired_alerts"][:8])
        print(f"{name:<16} {mult:5.3f} {n_fire:5d} {n_act:4d} {n_cs:5d} "
              f"{load:7.3f}  {aids}")

    # Detailed view for IAA and Quercetin
    print("\n--- IAA detailed alerts ---")
    mol = Chem.MolFromSmiles("OC(=O)Cc1c[nH]c2ccccc12")
    alerts = detect_alerts(mol)
    for a in alerts:
        spec = ALERT_LIBRARY[a["alert_id"]]
        pen = a["severity"] * a["probability"] * (1.0 - a["recovery"])
        tag = " [ACTIVATING]" if _is_activating(spec) else ""
        print(f"  {a['alert_id']:<10} pen={pen:.3f} dir={a['direction']:<14}{tag}  {a['description'][:60]}")

    print("\n--- Quercetin detailed alerts ---")
    mol = Chem.MolFromSmiles("O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12")
    alerts = detect_alerts(mol)
    mult, det = compute_gate5_multiplier(mol)
    for a in alerts:
        pen = det["penalties"].get(a["alert_id"], 0)
        tag = " [ACTIVATING]" if a["alert_id"] in det["activating_alerts"] else ""
        print(f"  {a['alert_id']:<10} pen={pen:.3f} dir={a['direction']:<14}{tag}  {a['description'][:60]}")
    print(f"  → multiplier = {mult:.3f}")

    print("\n--- TNT detailed alerts ---")
    mol = Chem.MolFromSmiles("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]")
    alerts = detect_alerts(mol)
    mult, det = compute_gate5_multiplier(mol)
    for a in alerts:
        pen = det["penalties"].get(a["alert_id"], 0)
        cs = " [CLASS_SWITCH]" if a["gate3_class_switch"] else ""
        print(f"  {a['alert_id']:<10} pen={pen:.3f} dir={a['direction']:<14}{cs}  {a['description'][:60]}")
    print(f"  → multiplier = {mult:.3f}")
