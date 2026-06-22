"""
Scoring gates for the Plant Systemic Bioavailability cascade.

Implements Gates 1-4 and the Composite Integrator from the Master
Checkpoint v2 (Section 3).  Each gate is a pure function:
descriptor dict + classification dict in, (score, details) out.

The classification dict is produced by the Module 7 Class Classifier
and must contain at minimum:
    module7_class   : str   — one of the 18 base classes
    gate3_mechanism : str   — routing key for Gate 3
    carrier_alpha   : float — carrier recognition probability [0,1]
"""

import math

# ── Helpers ────────────────────────────────────────────────────────────

def _g(x, mu, sigma):
    """Gaussian desirability centered at mu with width sigma."""
    return math.exp(-0.5 * ((x - mu) / sigma) ** 2)


def _dval(d, key):
    """Extract numeric value from a descriptor dict entry."""
    v = d.get(key)
    if v is None:
        return None
    if isinstance(v, dict):
        return v.get("value")
    return v


def _clip(x, lo=0.0, hi=1.0):
    return max(lo, min(hi, x))


# ── Gate 1: Cuticle Penetration ───────────────────────────────────────

# pp-LFER coefficients for isolated astomatous cuticles (Qi et al. 2020)
_PPLFER = {"c": 0.52, "e": 0.60, "s": -0.72,
           "a": -1.60, "b": -3.58, "v": 3.62}


def score_gate1(descriptors, classification=None, use_ppLFER=False):
    """Gate 1: Cuticle penetration score (0-1).

    Parameters
    ----------
    descriptors : dict
        Output of compute_descriptors (key → {value, provenance}).
    classification : dict, optional
        Unused by Gate 1; kept for interface consistency.
    use_ppLFER : bool
        If True and Abraham descriptors are available, use the
        pp-LFER model instead of the default Gaussian.

    Returns
    -------
    (float, dict) — (score, details)
    """
    logKow = _dval(descriptors, "logKow_value") or 0.0
    logD = _dval(descriptors, "logD_pH5p5")
    if logD is None:
        logD = logKow
    mw = _dval(descriptors, "mw") or 200.0
    sol = _dval(descriptors, "water_sol_mgL")
    mp = _dval(descriptors, "melting_point_C")
    vp = _dval(descriptors, "vapor_pressure_mPa")

    details = {"mode": "gaussian", "logD_used": round(logD, 2)}

    # pp-LFER mode
    if use_ppLFER:
        E = _dval(descriptors, "abraham_E")
        S = _dval(descriptors, "abraham_S")
        A = _dval(descriptors, "abraham_A")
        B = _dval(descriptors, "abraham_B")
        V = _dval(descriptors, "abraham_V")
        if all(v is not None for v in (E, S, A, B, V)):
            c = _PPLFER
            log_K = (c["c"] + c["e"] * E + c["s"] * S +
                     c["a"] * A + c["b"] * B + c["v"] * V)
            score = _clip(log_K / 6.0)
            score = max(score, 0.05) if logKow < -1 else score
            details.update(mode="ppLFER", log_Kcuticle=round(log_K, 3))
            return (round(score, 4), details)

    # Default Gaussian mode
    f_logKow = _g(logD, 2.5, 1.5)
    f_MW = _g(mw, 300, 150) if mw <= 800 else 0.05
    f_Sol = _clip(math.log10(sol + 1) / 3.0) if sol is not None else 0.5
    f_MP = _g(mp, 150, 100) if mp is not None else 0.5
    f_VP = (0.2 if vp > 10 else 1.0) if vp is not None else 1.0

    score = (0.45 * f_logKow + 0.25 * f_MW + 0.20 * f_Sol +
             0.05 * f_MP + 0.05 * f_VP)

    # Polar pathway floor
    if logKow < -1:
        score = max(score, 0.05)

    details.update(f_logKow=round(f_logKow, 4), f_MW=round(f_MW, 4),
                   f_Sol=round(f_Sol, 4), f_MP=round(f_MP, 4),
                   f_VP=round(f_VP, 4))
    return (round(_clip(score), 4), details)


# ── Gate 2: Phloem Loading ────────────────────────────────────────────

_CARRIER_SCORE_MAP = {
    "carrier_override":   0.85,
    "inverse_kleier":     0.20,
    "efflux_penalty":     0.10,
    "negative_control":   0.05,
    "full_kleier":        0.50,
    "partial_kleier":     0.60,
    "passive_moderate":   0.45,
    "intermediate_perm":  0.55,
    "retention_penalized": 0.35,
    "wax_bound":          0.05,
    "permanent_dication": 0.02,
}


def score_gate2(descriptors, classification):
    """Gate 2: Phloem loading score (0-1).

    Parameters
    ----------
    descriptors : dict
        Descriptor dict.
    classification : dict
        Must contain carrier_alpha and gate3_mechanism.

    Returns
    -------
    (float, dict) — (score, details)
    """
    alpha = classification.get("carrier_alpha", 0.0)
    mechanism = classification.get("gate3_mechanism", "full_kleier")

    pka_acid = _dval(descriptors, "pka_acidic")
    mw = _dval(descriptors, "mw") or 200.0
    tpsa = _dval(descriptors, "tpsa") or 60.0

    # Passive sub-score
    if pka_acid is not None:
        f_pKa = _g(pka_acid, 5.0, 1.5)
    else:
        f_pKa = 0.30  # neutrals / bases

    if 100 <= mw <= 400:
        f_MW = 1.0
    elif mw < 100:
        f_MW = _clip(mw / 100.0)
    else:
        f_MW = _clip(1.0 - (mw - 400) / 400.0, 0.05, 1.0)

    f_TPSA = _g(tpsa, 60, 40)

    passive = 0.40 * f_pKa + 0.30 * f_MW + 0.30 * f_TPSA

    # Carrier sub-score
    carrier = _CARRIER_SCORE_MAP.get(mechanism, 0.50)

    score = (1.0 - alpha) * passive + alpha * carrier

    details = {
        "passive_score": round(passive, 4),
        "carrier_score": round(carrier, 4),
        "carrier_alpha": round(alpha, 3),
        "f_pKa": round(f_pKa, 4),
        "f_MW": round(f_MW, 4),
        "f_TPSA": round(f_TPSA, 4),
    }
    return (round(_clip(score), 4), details)


# ── Gate 3: Vascular Transport (2-D Kleier Surface) ──────────────────

def _full_kleier_score(logKow, pKa):
    """Classic 2-D Kleier: norm(logCf) × R(logKow)."""
    logP_nitella = 1.20 * logKow - 5.86
    logCf = logP_nitella + 2.303 * (7.5 - pKa)
    norm_logCf = _clip((logCf + 4.0) / 6.0)
    R = math.exp(-((logKow - 1.5) / 1.5) ** 2)
    return norm_logCf * R, logCf, R


def score_gate3(descriptors, classification):
    """Gate 3: Vascular transport score (0-1).

    Mechanism-branched based on classification.gate3_mechanism.

    Returns
    -------
    (float, dict) — (score, details including mechanism, logCf, R_logKow)
    """
    mechanism = classification.get("gate3_mechanism", "full_kleier")
    alpha = classification.get("carrier_alpha", 0.0)
    logKow = _dval(descriptors, "logKow_value") or 0.0
    pka_acid = _dval(descriptors, "pka_acidic")
    pka_base = _dval(descriptors, "pka_basic")

    details = {"mechanism": mechanism, "logKow": round(logKow, 2)}
    score = 0.0

    if mechanism == "full_kleier":
        if pka_acid is None:
            pka_acid = 5.0  # assume mild acid if classified as weak acid
        raw, logCf, R = _full_kleier_score(logKow, pka_acid)
        score = raw
        details.update(logCf=round(logCf, 2), R_logKow=round(R, 4),
                       pKa_used=pka_acid)

    elif mechanism == "carrier_override":
        score = _clip(alpha, 0.50, 0.90)
        details["carrier_alpha"] = round(alpha, 3)

    elif mechanism == "inverse_kleier":
        pKa = pka_base if pka_base is not None else 8.0
        logP_n = 1.20 * logKow - 5.86
        logCf_phloem = logP_n + 2.303 * (pKa - 7.5)
        score = _clip(0.05 + 0.15 * max(0, -logCf_phloem / 5.0), 0, 0.20)
        details.update(logCf_phloem=round(logCf_phloem, 2), pKa_used=pKa)

    elif mechanism == "partial_kleier":
        pKa_eff = pka_acid if pka_acid is not None else 5.0
        kleier_raw, logCf, R = _full_kleier_score(logKow, pKa_eff)
        score = 0.5 * kleier_raw + 0.5 * _clip(alpha, 0.3, 0.8)
        details.update(logCf=round(logCf, 2), R_logKow=round(R, 4),
                       carrier_alpha=round(alpha, 3))

    elif mechanism == "passive_moderate":
        score = 0.3 + 0.3 * _g(logKow, 1.5, 2.0)
        details["logKow_desirability"] = round(_g(logKow, 1.5, 2.0), 4)

    elif mechanism == "negative_control":
        score = 0.02

    elif mechanism == "efflux_penalty":
        score = max(0.01, 0.15 + alpha)
        details["carrier_alpha"] = round(alpha, 3)

    elif mechanism == "intermediate_perm":
        logKow_d = _g(logKow, 0.0, 1.0)
        score = 0.3 + 0.3 * logKow_d + 0.1 * _clip(alpha)
        details.update(logKow_desirability=round(logKow_d, 4),
                       carrier_alpha=round(alpha, 3))

    elif mechanism == "retention_penalized":
        R = math.exp(-((logKow - 1.5) / 1.5) ** 2)
        score = 0.15 * R
        details["R_logKow"] = round(R, 4)

    elif mechanism == "wax_bound":
        score = 0.02

    elif mechanism == "permanent_dication":
        score = 0.01

    else:
        R = math.exp(-((logKow - 1.5) / 1.5) ** 2)
        score = 0.10 * R
        details["R_logKow"] = round(R, 4)
        details["fallback"] = True

    details["score"] = round(_clip(score), 4)
    return (round(_clip(score), 4), details)


# ── Gate 4: Phloem Unloading ─────────────────────────────────────────

def score_gate4(descriptors, classification, gate5_alerts=None):
    """Gate 4: Phloem unloading score (0-1).

    Parameters
    ----------
    descriptors : dict
        Descriptor dict.
    classification : dict
        Module 7 classification.
    gate5_alerts : list[str], optional
        Fired Gate 5 alert IDs.  B12 (malonyl-cap) strengthens
        vacuolar coupling and competes with unloading.

    Returns
    -------
    (float, dict) — (score, details)
    """
    mechanism = classification.get("gate3_mechanism", "full_kleier")

    # Short-circuit for clearly non-phloem classes
    if mechanism in ("negative_control", "efflux_penalty",
                     "wax_bound", "permanent_dication"):
        return (0.05, {"mechanism": mechanism, "short_circuit": True})

    # Carrier-mediated unloading assumed for override classes
    if mechanism == "carrier_override":
        return (0.65, {"mechanism": mechanism, "carrier_unloading": True})

    pka_acid = _dval(descriptors, "pka_acidic")
    pka_base = _dval(descriptors, "pka_basic")
    tpsa = _dval(descriptors, "tpsa") or 60.0
    hbd = _dval(descriptors, "hbd") or 0
    mw = _dval(descriptors, "mw") or 200.0

    # f_pKa: acid release at sink apoplast pH 4.5-5.5
    if pka_acid is not None:
        f_pKa = _g(pka_acid, 4.5, 1.5)
    elif pka_base is not None:
        f_pKa = 0.20  # bases are vacuole-trapped
    else:
        f_pKa = 0.40  # neutrals

    # f_polar: polar recognition
    f_polar = _clip(tpsa / 60.0 + hbd * 0.15)

    # f_MW
    if 100 <= mw <= 350:
        f_MW = 1.0
    elif mw < 100:
        f_MW = _clip(mw / 100.0)
    else:
        f_MW = _clip(1.0 - (mw - 350) / 350.0, 0.1, 1.0)

    # f_stability
    f_stability = 1.0
    if gate5_alerts and "B12" in gate5_alerts:
        # Malonyl-cap strengthens vacuolar sequestration
        f_stability = 0.80

    score = 0.40 * f_pKa + 0.25 * f_polar + 0.20 * f_MW + 0.15 * f_stability

    details = {
        "f_pKa": round(f_pKa, 4),
        "f_polar": round(f_polar, 4),
        "f_MW": round(f_MW, 4),
        "f_stability": round(f_stability, 4),
        "malonyl_cap": gate5_alerts is not None and "B12" in (gate5_alerts or []),
    }
    return (round(_clip(score), 4), details)


# ── Composite Integrator ─────────────────────────────────────────────

def score_composite(gate_scores, gate5_multiplier=1.0,
                    descriptors=None, classification=None):
    """Composite = G1 × G2 × G3 × G4 × G5_multiplier.

    Parameters
    ----------
    gate_scores : dict
        {"gate1": float, "gate2": float, "gate3": float, "gate4": float}
    gate5_multiplier : float
        Product of (1 - penalty_i) from fired Gate 5 alerts.
    descriptors : dict, optional
        For diagnostic flag generation.
    classification : dict, optional
        For diagnostic flag generation.

    Returns
    -------
    dict with composite, gate_scores, bottleneck, design_flags, etc.
    """
    g1 = gate_scores.get("gate1", 0.5)
    g2 = gate_scores.get("gate2", 0.5)
    g3 = gate_scores.get("gate3", 0.5)
    g4 = gate_scores.get("gate4", 0.5)
    g5 = _clip(gate5_multiplier)

    composite = g1 * g2 * g3 * g4 * g5

    # Bottleneck
    named = {"gate1": g1, "gate2": g2, "gate3": g3,
             "gate4": g4, "gate5": g5}
    bottleneck = min(named, key=named.get)

    # Design flags
    flags = []

    logKow = None
    pka_acid = None

    if descriptors:
        logKow = _dval(descriptors, "logKow_value")
        pka_acid = _dval(descriptors, "pka_acidic")

    if g1 < 0.15:
        flags.append("GATE1_CUTICLE_BLOCK")

    if logKow is not None and logKow < -2:
        flags.append("GATE1_ADJUVANT_DEPENDENT")

    if g5 < 0.50:
        flags.append("GATE5_METABOLIC_RISK")

    if logKow is not None and logKow > 3.5 and g3 < 0.30:
        flags.append("GATE3_LEAKY_LIPOPHILE")

    if pka_acid is not None and 3.0 <= pka_acid <= 6.0 and g3 > 0.50:
        flags.append("AMBIMOBILE")

    if classification:
        cls = classification.get("module7_class", "")
        if cls == "TERTIARY_AMINE_ALKALOID":
            flags.append("INVERSE_KLEIER")
        if cls in ("STEROID_HORMONE", "HIGH_LIPOPHILE"):
            flags.append("NON_SYSTEMIC")
        if cls == "OTHER_permanent_dication":
            flags.append("PERMANENT_DICATION")

        if classification.get("has_chiral_centers", False):
            flags.append("STEREO_UNRESOLVED")

    return {
        "composite": round(composite, 4),
        "gate_scores": {k: round(v, 4) for k, v in named.items()},
        "bottleneck": bottleneck,
        "bottleneck_score": round(named[bottleneck], 4),
        "design_flags": flags,
    }


# ── Full scoring convenience ─────────────────────────────────────────

def score_all_gates(descriptors, classification,
                    gate5_alerts=None, gate5_multiplier=1.0,
                    use_ppLFER=False):
    """Run the full cascade and return composite + per-gate details.

    Returns
    -------
    dict with:
        composite_result : dict from score_composite
        gate1 : (score, details)
        gate2 : (score, details)
        gate3 : (score, details)
        gate4 : (score, details)
        gate5_multiplier : float
    """
    g1_score, g1_det = score_gate1(descriptors, classification, use_ppLFER)
    g2_score, g2_det = score_gate2(descriptors, classification)
    g3_score, g3_det = score_gate3(descriptors, classification)
    g4_score, g4_det = score_gate4(descriptors, classification, gate5_alerts)

    gate_scores = {
        "gate1": g1_score,
        "gate2": g2_score,
        "gate3": g3_score,
        "gate4": g4_score,
    }

    comp = score_composite(gate_scores, gate5_multiplier,
                           descriptors, classification)

    return {
        "composite_result": comp,
        "gate1": (g1_score, g1_det),
        "gate2": (g2_score, g2_det),
        "gate3": (g3_score, g3_det),
        "gate4": (g4_score, g4_det),
        "gate5_multiplier": gate5_multiplier,
    }


# ── Demo / test ──────────────────────────────────────────────────────

if __name__ == "__main__":
    from rdkit import Chem
    from plant_scoring.descriptors import compute_descriptors

    # Seed compounds with expected Module 7 classifications
    SEEDS = [
        ("OC(=O)CNCP(=O)(O)O", "Glyphosate",
         {"module7_class": "PHOSPHONATE_NUTRIENT_MIMIC",
          "gate3_mechanism": "carrier_override", "carrier_alpha": 0.80,
          "has_chiral_centers": False},
         "3_phloem"),
        ("OC(=O)COc1ccc(Cl)cc1Cl", "2,4-D",
         {"module7_class": "WEAK_ACID_ORGANIC",
          "gate3_mechanism": "full_kleier", "carrier_alpha": 0.0,
          "has_chiral_centers": False},
         "3_phloem"),
        ("OC(=O)Cc1c[nH]c2ccccc12", "IAA",
         {"module7_class": "WEAK_ACID_ORGANIC",
          "gate3_mechanism": "full_kleier", "carrier_alpha": 0.0,
          "has_chiral_centers": False},
         "3_phloem"),
        ("O=[N+]([O-])NC1=NCCN1Cc1ccc(Cl)nc1", "Imidacloprid",
         {"module7_class": "POLAR_NEUTRAL_SYSTEMIC",
          "gate3_mechanism": "intermediate_perm", "carrier_alpha": 0.30,
          "has_chiral_centers": False},
         "3_phloem"),
        ("CC(C)(C)C(O)(Cn1cncn1)CCc1ccc(Cl)cc1", "Tebuconazole",
         {"module7_class": "NEUTRAL_LIPOPHILE",
          "gate3_mechanism": "retention_penalized", "carrier_alpha": 0.0,
          "has_chiral_centers": True},
         "2_xylem"),
        ("CN1CCC[C@H]1c1cccnc1", "Nicotine",
         {"module7_class": "TERTIARY_AMINE_ALKALOID",
          "gate3_mechanism": "inverse_kleier", "carrier_alpha": 0.0,
          "has_chiral_centers": True},
         "2_xylem"),
        ("C[n+]1ccc(-c2cc[n+](C)cc2)cc1", "Paraquat",
         {"module7_class": "OTHER_permanent_dication",
          "gate3_mechanism": "permanent_dication", "carrier_alpha": 0.0,
          "has_chiral_centers": False},
         "0_contact"),
        ("N#Cc1c(Cl)c(C#N)c(Cl)c(Cl)c1Cl", "Chlorothalonil",
         {"module7_class": "NEUTRAL_LIPOPHILE",
          "gate3_mechanism": "retention_penalized", "carrier_alpha": 0.0,
          "has_chiral_centers": False},
         "0_contact"),
    ]

    print("=" * 115)
    print("PLANT SYSTEMIC BIOAVAILABILITY — GATE SCORING DEMO")
    print("=" * 115)
    hdr = (f"{'Name':<16} {'Class':<14} {'G1':>6} {'G2':>6} {'G3':>6} "
           f"{'G4':>6} {'G5':>5} {'Comp':>6} {'Neck':<8} {'Expected':<10} "
           f"{'Flags'}")
    print(hdr)
    print("-" * 115)

    for smi, name, cls, expected in SEEDS:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"{name:<16} INVALID SMILES")
            continue
        desc = compute_descriptors(mol)
        r = score_all_gates(desc, cls, gate5_multiplier=1.0)

        c = r["composite_result"]
        flags = ", ".join(c["design_flags"]) if c["design_flags"] else "---"
        print(f"{name:<16} {cls['module7_class'][:13]:<14} "
              f"{r['gate1'][0]:6.3f} {r['gate2'][0]:6.3f} "
              f"{r['gate3'][0]:6.3f} {r['gate4'][0]:6.3f} "
              f"{r['gate5_multiplier']:5.2f} {c['composite']:6.4f} "
              f"{c['bottleneck']:<8} {expected:<10} {flags}")

    print()
    print("Scoring interpretation:")
    print("  3_phloem compounds should have composite > 0.05")
    print("  2_xylem compounds should have low Gate 3 (inverse Kleier / retention penalty)")
    print("  0_contact compounds should have composite ≈ 0 (multiple gate failures)")
