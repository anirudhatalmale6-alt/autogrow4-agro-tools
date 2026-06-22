#!/usr/bin/env python3
"""
Structure standardization and identifier resolution for the Plant Systemic
Bioavailability scoring pipeline (Phase B1 / R4).

Provides:
  standardize_smiles()   — canonical SMILES, salt-strip, tautomer, InChIKey
  resolve_protonation()  — net charge & ionization state at a given pH
  batch_standardize()    — process a list of SMILES
  fetch_identifiers()    — PubChem PUG-REST lookup (CID, CAS, DTXSID, …)
"""

import re
import time
from typing import Optional

import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
from rdkit.Chem.MolStandardize import rdMolStandardize

# ---------------------------------------------------------------------------
# Ionizable-group library (SMARTS → typical pKa)
# ---------------------------------------------------------------------------

ACIDIC_GROUPS = [
    ("sulfonic_acid",   "[SX4](=O)(=O)[OH]",            -1.0),
    ("phosphonic_1",    "[PX4](=O)([OH])([OH])",          2.2),
    ("carboxylic_acid", "[CX3](=O)[OH]",                  4.0),
    ("tetrazole",       "c1nnn[nH]1",                     4.7),
    ("phosphonic_2",    "[PX4](=O)([OH])([O-,OH])",       5.5),
    ("thiol",           "[SX2H]",                         8.5),
    ("sulfonamide",     "[NX3;H1][SX4](=O)=O",          10.0),
    ("phenol",          "[cR1][OH]",                     10.0),
]

BASIC_GROUPS = [
    # 1,2,4-triazole: weak base, pKa ~2.8 (N-substituted or NH)
    ("triazole_124",    "[n]1[c][n][c][n]1",              2.8),
    # Pyridine: exclude quaternary (already permanently charged)
    ("pyridine",        "[nX2;!$([n+])]1ccccc1",          5.2),
    ("imidazole",       "c1cnc[nH]1",                     6.0),
    ("tertiary_amine",  "[NX3;H0;!$(N-[!#6]);!$(N-C=O)]([CX4])([CX4])[CX4]", 9.8),
    ("primary_amine",   "[NX3;H2;!$(N-C=O);!$(N~[!#6])]",10.5),
    ("secondary_amine", "[NX3;H1;!$(N-C=O);!$(N~[!#6])]([CX4])[CX4]", 10.7),
    # Guanidine: exclude nitroguanidines (imidacloprid) where N is bonded
    # to [N+](=O) — the nitro drops pKa to ~1, effectively non-basic
    ("guanidine",       "[NX3;!$(N-[N+])][CX3](=[NX2])[NX3;!$(N-[N+])]",
                                                         12.5),
]

# Pre-compile SMARTS once at import
_ACIDIC_PATS = [(n, Chem.MolFromSmarts(s), pka) for n, s, pka in ACIDIC_GROUPS]
_BASIC_PATS  = [(n, Chem.MolFromSmarts(s), pka) for n, s, pka in BASIC_GROUPS]

# Reusable standardization singletons
_LFC      = rdMolStandardize.LargestFragmentChooser()
_UNCHARGER = rdMolStandardize.Uncharger()
_TAUT_ENUM = rdMolStandardize.TautomerEnumerator()

# PubChem rate-limit state
_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_PUBCHEM_LAST_CALL = 0.0
_PUBCHEM_MIN_INTERVAL = 0.25  # 4 req/s max; PUG-REST allows 5/s


# ---------------------------------------------------------------------------
# Core: standardize_smiles
# ---------------------------------------------------------------------------

def standardize_smiles(smiles: str) -> dict:
    """Standardize a SMILES string through the full pipeline.

    Returns a dict with keys: canonical_smiles, original_smiles, inchikey,
    mol, standardization_log, is_valid.
    """
    result = {
        "original_smiles": smiles,
        "canonical_smiles": None,
        "inchikey": None,
        "mol": None,
        "standardization_log": [],
        "is_valid": False,
    }
    log = result["standardization_log"]

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        log.append("FAILED: invalid SMILES")
        return result

    original_smi = Chem.MolToSmiles(mol)

    # 1. Cleanup (valence, charge normalization)
    mol = rdMolStandardize.Cleanup(mol)
    log.append("cleanup")

    # 2. Normalize
    mol = rdMolStandardize.Normalize(mol)
    log.append("normalize")

    # 3. Largest fragment (salt strip)
    num_frags = len(Chem.GetMolFrags(mol))
    mol = _LFC.choose(mol)
    if num_frags > 1:
        log.append(f"salt_stripped ({num_frags} frags → 1)")
    else:
        log.append("single_fragment")

    # 4. Neutralize
    pre_neutralize = Chem.MolToSmiles(mol)
    mol = _UNCHARGER.uncharge(mol)
    if Chem.MolToSmiles(mol) != pre_neutralize:
        log.append("neutralized")

    # 5. Tautomer canonicalization
    pre_taut = Chem.MolToSmiles(mol)
    mol = _TAUT_ENUM.Canonicalize(mol)
    if Chem.MolToSmiles(mol) != pre_taut:
        log.append("tautomer_resolved")

    canonical = Chem.MolToSmiles(mol)

    # InChIKey
    inchi = MolToInchi(mol)
    inchikey = InchiToInchiKey(inchi) if inchi else None

    result["canonical_smiles"] = canonical
    result["inchikey"] = inchikey
    result["mol"] = mol
    result["is_valid"] = True
    return result


# ---------------------------------------------------------------------------
# Protonation state estimation
# ---------------------------------------------------------------------------

def detect_ionizable_groups(mol) -> list[dict]:
    """Find all ionizable functional groups and their estimated pKa values."""
    groups = []

    for name, pat, pka in _ACIDIC_PATS:
        if pat is None:
            continue
        matches = mol.GetSubstructMatches(pat)
        for i, match in enumerate(matches):
            groups.append({
                "group": name,
                "type": "acidic",
                "pKa": pka,
                "atoms": list(match),
                "instance": i,
            })

    for name, pat, pka in _BASIC_PATS:
        if pat is None:
            continue
        matches = mol.GetSubstructMatches(pat)
        for i, match in enumerate(matches):
            groups.append({
                "group": name,
                "type": "basic",
                "pKa": pka,
                "atoms": list(match),
                "instance": i,
            })

    return groups


def resolve_protonation(mol, pH: float) -> dict:
    """Estimate net charge and dominant ionization state at a given pH.

    Henderson-Hasselbalch:
      acid deprotonated fraction: α = 10^(pH-pKa) / (1 + 10^(pH-pKa))
      base protonated fraction:   α = 10^(pKa-pH) / (1 + 10^(pKa-pH))

    Returns dict with: net_charge, groups (with fraction_ionized), pH.
    """
    groups = detect_ionizable_groups(mol)
    net_charge = 0
    annotated = []

    for g in groups:
        pka = g["pKa"]
        if g["type"] == "acidic":
            # fraction deprotonated (bearing negative charge)
            exp = min(pH - pka, 50)
            alpha = 10**exp / (1 + 10**exp)
            charge_contribution = -alpha
        else:
            # fraction protonated (bearing positive charge)
            exp = min(pka - pH, 50)
            alpha = 10**exp / (1 + 10**exp)
            charge_contribution = +alpha

        net_charge += charge_contribution
        annotated.append({
            **g,
            "fraction_ionized": round(alpha, 4),
            "charge_contribution": round(charge_contribution, 4),
        })

    # Permanent charges: quaternary sp3 N+, aromatic n+ (paraquat/diquat),
    # permanent anions (sulfonate, etc.)
    _perm_pos_pats = [
        Chem.MolFromSmarts("[NX4;+1;!$([NH4+])]"),     # quaternary N
        Chem.MolFromSmarts("[n+]"),                      # aromatic N+ (paraquat)
    ]
    _perm_neg_pats = [
        Chem.MolFromSmarts("[OX1;-1]~[#16X4]"),         # sulfonate
        Chem.MolFromSmarts("[OX1;-1]~[#15X4]"),         # phosphonate anion
    ]
    permanent_pos = sum(len(mol.GetSubstructMatches(p))
                        for p in _perm_pos_pats if p)
    permanent_neg = sum(len(mol.GetSubstructMatches(p))
                        for p in _perm_neg_pats if p)
    net_charge += permanent_pos - permanent_neg

    return {
        "pH": pH,
        "net_charge": round(net_charge, 3),
        "net_charge_int": round(net_charge),
        "groups": annotated,
        "permanent_positive": permanent_pos,
        "permanent_negative": permanent_neg,
    }


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------

def batch_standardize(smiles_list: list[str]) -> list[dict]:
    """Standardize a list of SMILES strings. Each result dict also gets
    protonation data at pH 5.5 (apoplast) and 7.5 (sieve tube)."""
    results = []
    for smi in smiles_list:
        rec = standardize_smiles(smi)
        if rec["is_valid"]:
            mol = rec["mol"]
            rec["protonation_pH5_5"] = resolve_protonation(mol, 5.5)
            rec["protonation_pH7_5"] = resolve_protonation(mol, 7.5)
        else:
            rec["protonation_pH5_5"] = None
            rec["protonation_pH7_5"] = None
        # Drop the Mol object for serialization-friendliness
        rec["mol_present"] = rec["mol"] is not None
        results.append(rec)
    return results


# ---------------------------------------------------------------------------
# PubChem identifier resolution
# ---------------------------------------------------------------------------

def _pubchem_throttle():
    global _PUBCHEM_LAST_CALL
    elapsed = time.time() - _PUBCHEM_LAST_CALL
    if elapsed < _PUBCHEM_MIN_INTERVAL:
        time.sleep(_PUBCHEM_MIN_INTERVAL - elapsed)
    _PUBCHEM_LAST_CALL = time.time()


def _pubchem_get(url: str, timeout: int = 15) -> Optional[dict]:
    _pubchem_throttle()
    try:
        r = requests.get(url, timeout=timeout)
        if r.status_code == 200:
            return r.json()
    except (requests.RequestException, ValueError):
        pass
    return None


_CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
_DTXSID_RE = re.compile(r"^DTXSID\d+$")


def fetch_identifiers(query: str, use_pubchem: bool = True) -> Optional[dict]:
    """Look up a compound by name or SMILES via PubChem PUG-REST.

    Returns dict with: cid, canonical_smiles, inchikey, common_name,
    cas_rn, dtxsid.  Returns None if PubChem is unreachable or the
    compound is not found.
    """
    if not use_pubchem:
        return None

    # Determine query type: if it parses as a molecule, search by SMILES
    mol = Chem.MolFromSmiles(query)
    if mol is not None:
        search_type = "smiles"
        search_val = Chem.MolToSmiles(mol)
    else:
        search_type = "name"
        search_val = query

    # Step 1: get properties
    prop_url = (
        f"{_PUBCHEM_BASE}/compound/{search_type}/{requests.utils.quote(search_val)}"
        f"/property/CanonicalSMILES,InChIKey,MolecularFormula,IUPACName/JSON"
    )
    prop_data = _pubchem_get(prop_url)
    if prop_data is None:
        return None

    try:
        props = prop_data["PropertyTable"]["Properties"][0]
    except (KeyError, IndexError):
        return None

    cid = props.get("CID")
    result = {
        "cid": cid,
        "canonical_smiles": props.get("CanonicalSMILES") or props.get("ConnectivitySMILES"),
        "inchikey": props.get("InChIKey"),
        "iupac_name": props.get("IUPACName"),
        "molecular_formula": props.get("MolecularFormula"),
        "common_name": None,
        "cas_rn": None,
        "dtxsid": None,
    }

    if cid is None:
        return result

    # Step 2: get synonyms for common name, CAS, DTXSID
    syn_url = f"{_PUBCHEM_BASE}/compound/cid/{cid}/synonyms/JSON"
    syn_data = _pubchem_get(syn_url)
    if syn_data is not None:
        try:
            syns = syn_data["InformationList"]["Information"][0]["Synonym"]
        except (KeyError, IndexError):
            syns = []

        if syns:
            result["common_name"] = syns[0]

        for s in syns:
            if result["cas_rn"] is None and _CAS_RE.match(s):
                result["cas_rn"] = s
            if result["dtxsid"] is None and _DTXSID_RE.match(s):
                result["dtxsid"] = s
            if result["cas_rn"] and result["dtxsid"]:
                break

    return result


# ---------------------------------------------------------------------------
# Demo / self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    SEED_COMPOUNDS = {
        "Glyphosate":      "OC(=O)CNCP(=O)(O)O",
        "2,4-D":           "OC(=O)COc1ccc(Cl)cc1Cl",
        "IAA":             "OC(=O)Cc1c[nH]c2ccccc12",
        "Imidacloprid":    "O=[N+]([O-])/N=C1\\NCCN1Cc1ccc(Cl)nc1",
        "Tebuconazole":    "CC(C)(C)C(O)(Cn1cncn1)CCc1ccc(Cl)cc1",
        "Nicotine":        "CN1CCC[C@H]1c1cccnc1",
        "Paraquat":        "C[n+]1ccc(-c2cc[n+](C)cc2)cc1",
        "Chlorothalonil":  "N#Cc1c(Cl)c(C#N)c(Cl)c(Cl)c1Cl",
    }

    print("=" * 90)
    print("STRUCTURE STANDARDIZATION MODULE — SELF-TEST")
    print("=" * 90)

    for name, smi in SEED_COMPOUNDS.items():
        print(f"\n{'─' * 90}")
        print(f"  {name}  ({smi})")
        print(f"{'─' * 90}")

        rec = standardize_smiles(smi)
        print(f"  Valid:      {rec['is_valid']}")
        print(f"  Canonical:  {rec['canonical_smiles']}")
        print(f"  InChIKey:   {rec['inchikey']}")
        print(f"  Log:        {', '.join(rec['standardization_log'])}")

        if rec["is_valid"]:
            mol = rec["mol"]
            for pH_label, pH_val in [("apoplast pH 5.5", 5.5),
                                      ("sieve tube pH 7.5", 7.5)]:
                prot = resolve_protonation(mol, pH_val)
                print(f"  {pH_label}: net_charge = {prot['net_charge']:+.2f} "
                      f"(rounded: {prot['net_charge_int']:+d})")
                for g in prot["groups"]:
                    print(f"    {g['group']:20s} ({g['type']:6s}) "
                          f"pKa={g['pKa']:.1f}  "
                          f"α={g['fraction_ionized']:.3f}  "
                          f"Δq={g['charge_contribution']:+.3f}")
                if prot["permanent_positive"] or prot["permanent_negative"]:
                    print(f"    permanent: +{prot['permanent_positive']} "
                          f"/ -{prot['permanent_negative']}")

    # PubChem lookup demo (only first 3 to avoid hammering the API)
    print(f"\n{'=' * 90}")
    print("PUBCHEM IDENTIFIER LOOKUP")
    print(f"{'=' * 90}")

    for name in list(SEED_COMPOUNDS.keys())[:3]:
        print(f"\n  Querying: {name}")
        ids = fetch_identifiers(name)
        if ids:
            print(f"    CID:     {ids['cid']}")
            print(f"    CAS:     {ids['cas_rn']}")
            print(f"    DTXSID:  {ids['dtxsid']}")
            print(f"    Name:    {ids['common_name']}")
            print(f"    SMILES:  {ids['canonical_smiles']}")
        else:
            print(f"    (PubChem unreachable or not found)")
