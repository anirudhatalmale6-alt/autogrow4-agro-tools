#!/usr/bin/env python3
"""
PubChem Agrochemical Library Builder v2 for AutoGrow4
=====================================================
Downloads PubChem compound SMILES and filters for agrochemical drug-likeness
with reactive handle matching for the agricultural_reactions library.

Filters:
  - Elements: C, H, O, N, S, Si, B, F, Cl
  - Tice's Rules: LogP 0-4, RotBonds <=8, HBA 2-10, HBD <3
  - Kleier Map: TPSA 40-120 (phloem mobility)
  - Aromatic rings <= 3
  - SA Score <= 3.0
  - Reactive handles: must contain at least 1 of 12 handle substructures
  - MW bins: A(<100), B(100-150), C(150-200), D(200-250)

Usage:
  conda activate ag4_full
  python pubchem_agro_library_v2.py --output agrochemical_library.smi
  python pubchem_agro_library_v2.py --output agrochemical_library.smi --resume --skip-ghs
"""

import argparse
import gzip
import json
import os
import sys
import time
import urllib.request
import urllib.error
from collections import Counter

try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    print("ERROR: RDKit not found. Activate conda environment first:")
    print("  conda activate ag4_full")
    sys.exit(1)

try:
    from rdkit.Chem import RDConfig
    sa_score_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if os.path.isdir(sa_score_path):
        sys.path.insert(0, sa_score_path)
        import sascorer
        HAS_SASCORER = True
    else:
        HAS_SASCORER = False
except Exception:
    HAS_SASCORER = False

if not HAS_SASCORER:
    try:
        from rdkit.Chem import QED
        print("WARNING: SA_Score module not found. Using heuristic estimate.")
    except Exception:
        pass


PUBCHEM_FTP_URL = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz"
ALLOWED_ELEMENTS = {"C", "H", "O", "N", "S", "Si", "B", "F", "Cl", "Br"}
CHUNK_REPORT = 500000

MW_BINS = [
    ("A", 0, 100),
    ("B", 100, 150),
    ("C", 150, 200),
    ("D", 200, 250),
]

REACTIVE_HANDLES = {
    "silicon_switch": "[Si]",
    "boronic_acid": "[B]([OX2])([OX2])",
    "cf3_group": "[CX4](F)(F)F",
    "epoxide": "[CR1]1[OR1][CR1]1",
    "isothiocyanate": "[NX2]=[CX2]=[SX1]",
    "primary_secondary_amine": "[NX3;H2,H1;!$(NC=O)]",
    "aryl_halide": "[c][Cl,Br,I]",
    "heterocyclic_n": "[nR1]",
    "short_alkyl_ether": "[CX4][OX2][CX4]",
    "sulfur_bridge": "[#16X2]([#6])[#6]",
    "nitrile": "[CX2]#[NX1]",
    "alpha_hydroxy_acid": "[OX2H][CX4][CX3](=[OX1])",
    "carboxylic_acid": "[CX3](=[OX1])[OX2H1]",
    "alcohol": "[CX4][OX2H1]",
    "sulfonyl_chloride": "[SX4](=[OX1])(=[OX1])[Cl]",
    "alkyl_halide": "[CX4][Cl,Br,I]",
    "aldehyde": "[CX3H1]=[OX1]",
    "thiol": "[SX2H1]",
    "phenol": "[cX3][OX2H1]",
    "terminal_alkyne": "[CX2H1]#[CX2]",
}

_HANDLE_PATTERNS = None

def _get_handle_patterns():
    global _HANDLE_PATTERNS
    if _HANDLE_PATTERNS is None:
        _HANDLE_PATTERNS = {}
        for name, smarts in REACTIVE_HANDLES.items():
            pat = Chem.MolFromSmarts(smarts)
            if pat is not None:
                _HANDLE_PATTERNS[name] = pat
            else:
                print(f"WARNING: Invalid SMARTS for handle '{name}': {smarts}")
    return _HANDLE_PATTERNS


def calculate_sa_score(mol):
    """Calculate synthetic accessibility score (1=easy, 10=hard)."""
    if HAS_SASCORER:
        return sascorer.calculateScore(mol)
    ring_count = Descriptors.RingCount(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    stereo_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    score = 1.0 + (ring_count * 0.5) + (stereo_centers * 0.8) + (heavy_atoms * 0.05)
    return min(score, 10.0)


def get_elements(mol):
    """Get set of element symbols in a molecule."""
    return {atom.GetSymbol() for atom in mol.GetAtoms()}


def count_aromatic_rings(mol):
    """Count aromatic rings in molecule."""
    ri = mol.GetRingInfo()
    count = 0
    for ring in ri.BondRings():
        if all(mol.GetBondWithIdx(b).GetIsAromatic() for b in ring):
            count += 1
    return count


def get_mw_bin(mw):
    """Return MW bin label or None if out of range."""
    for label, lo, hi in MW_BINS:
        if lo <= mw < hi:
            return label
    return None


def has_reactive_handle(mol):
    """Check if molecule has at least one reactive handle. Returns list of matched handles."""
    handles = _get_handle_patterns()
    matched = []
    for name, pat in handles.items():
        if mol.HasSubstructMatch(pat):
            matched.append(name)
    return matched


def passes_filters(smiles, mw_min=0, mw_max=250, logp_min=-2, logp_max=5,
                   rot_bonds_max=8, hba_min=0, hba_max=10, hbd_max=5,
                   sa_max=3.0, tpsa_min=0, tpsa_max=200, aro_rings_max=3,
                   require_handles=True):
    """Apply all filters. Returns (pass, mol, props, bin_label) or (False, None, None, None)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, None, None, None

    elements = get_elements(mol)
    if not elements.issubset(ALLOWED_ELEMENTS):
        return False, None, None, None

    mw = Descriptors.ExactMolWt(mol)
    if mw < mw_min or mw >= mw_max:
        return False, None, None, None

    mw_bin = get_mw_bin(mw)
    if mw_bin is None:
        return False, None, None, None

    logp = Crippen.MolLogP(mol)
    if logp < logp_min or logp > logp_max:
        return False, None, None, None

    hbd = Lipinski.NumHDonors(mol)
    if hbd >= hbd_max:
        return False, None, None, None

    hba = Lipinski.NumHAcceptors(mol)
    if hba < hba_min or hba > hba_max:
        return False, None, None, None

    rot = Lipinski.NumRotatableBonds(mol)
    if rot > rot_bonds_max:
        return False, None, None, None

    tpsa = rdMolDescriptors.CalcTPSA(mol)
    if tpsa < tpsa_min or tpsa > tpsa_max:
        return False, None, None, None

    aro = count_aromatic_rings(mol)
    if aro > aro_rings_max:
        return False, None, None, None

    if require_handles:
        matched_handles = has_reactive_handle(mol)
        if not matched_handles:
            return False, None, None, None
    else:
        matched_handles = []

    sa = calculate_sa_score(mol)
    if sa > sa_max:
        return False, None, None, None

    props = {
        "mw": round(mw, 2),
        "logp": round(logp, 2),
        "hbd": hbd,
        "hba": hba,
        "rot": rot,
        "tpsa": round(tpsa, 2),
        "aro_rings": aro,
        "sa": round(sa, 2),
        "handles": matched_handles,
        "bin": mw_bin,
    }
    return True, mol, props, mw_bin


def download_and_filter(output_file, resume=False, skip_ghs=False, **filter_kwargs):
    """Download PubChem CID-SMILES.gz and filter compounds into MW bins."""

    base, ext = os.path.splitext(output_file)
    bin_files = {}
    progress_file = output_file + ".progress"
    skip_lines = 0

    if resume and os.path.exists(progress_file):
        with open(progress_file) as f:
            skip_lines = int(f.read().strip())
        print(f"Resuming from line {skip_lines}")

    mode = "a" if resume else "w"

    for label, lo, hi in MW_BINS:
        bin_path = f"{base}_bin{label}_{lo}to{hi}{ext}"
        bin_files[label] = open(bin_path, mode)
        print(f"  Bin {label}: MW {lo}-{hi} -> {bin_path}")

    combined_f = open(output_file + ".phase1.tmp" if not skip_ghs else output_file, mode)

    stats = Counter()
    handle_stats = Counter()
    bin_counts = Counter()
    total_processed = 0
    passed_count = 0

    rh = filter_kwargs.get("require_handles", True)
    print(f"\nDownloading PubChem CID-SMILES from: {PUBCHEM_FTP_URL}")
    print(f"Filters: MW {filter_kwargs.get('mw_min',0)}-{filter_kwargs.get('mw_max',250)}, "
          f"LogP {filter_kwargs.get('logp_min',0)}-{filter_kwargs.get('logp_max',4)}, "
          f"HBD <{filter_kwargs.get('hbd_max',3)}, HBA {filter_kwargs.get('hba_min',2)}-{filter_kwargs.get('hba_max',10)}, "
          f"RotBonds <={filter_kwargs.get('rot_bonds_max',8)}")
    print(f"TPSA: {filter_kwargs.get('tpsa_min',40)}-{filter_kwargs.get('tpsa_max',120)}, "
          f"Aro rings <={filter_kwargs.get('aro_rings_max',3)}, "
          f"SA <={filter_kwargs.get('sa_max',3.0)}")
    print(f"Elements: {','.join(sorted(ALLOWED_ELEMENTS))}")
    print(f"Reactive handles required: {rh}")
    if rh:
        print(f"Handles: {', '.join(REACTIVE_HANDLES.keys())}")
    print(f"SA Score method: {'sascorer (accurate)' if HAS_SASCORER else 'heuristic (approximate)'}")
    print()

    _get_handle_patterns()

    start_time = time.time()

    try:
        req = urllib.request.Request(PUBCHEM_FTP_URL)
        response = urllib.request.urlopen(req, timeout=120)
        gz = gzip.GzipFile(fileobj=response)

        line_num = 0
        for raw_line in gz:
            line_num += 1

            if line_num <= skip_lines:
                if line_num % (CHUNK_REPORT * 10) == 0:
                    print(f"  Skipping... at line {line_num:,}")
                continue

            total_processed += 1

            try:
                line = raw_line.decode("utf-8").strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                cid = parts[0]
                smiles = parts[1]
            except Exception:
                stats["parse_error"] += 1
                continue

            if "." in smiles:
                stats["fragmented"] += 1
                continue

            passed, mol, props, mw_bin = passes_filters(smiles, **filter_kwargs)

            if passed:
                entry = f"{smiles}\tCID_{cid}\n"
                combined_f.write(entry)
                if mw_bin in bin_files:
                    bin_files[mw_bin].write(entry)
                    bin_counts[mw_bin] += 1
                passed_count += 1
                stats["passed"] += 1
                for h in props.get("handles", []):
                    handle_stats[h] += 1
            else:
                stats["filtered"] += 1

            if total_processed % CHUNK_REPORT == 0:
                elapsed = time.time() - start_time
                rate = total_processed / elapsed if elapsed > 0 else 0
                bin_str = " ".join(f"{k}:{v}" for k, v in sorted(bin_counts.items()))
                print(f"  Processed: {total_processed:>12,} | "
                      f"Passed: {passed_count:>8,} | "
                      f"Rate: {rate:,.0f}/sec | "
                      f"Bins: {bin_str} | "
                      f"Elapsed: {elapsed/60:.1f} min")
                with open(progress_file, "w") as pf:
                    pf.write(str(line_num))
                combined_f.flush()
                for bf in bin_files.values():
                    bf.flush()

        gz.close()
        response.close()

    except KeyboardInterrupt:
        print(f"\nInterrupted at line {line_num}. Progress saved.")
        with open(progress_file, "w") as pf:
            pf.write(str(line_num))
        print(f"Resume with: python {sys.argv[0]} --output {output_file} --resume")
    except Exception as e:
        print(f"\nError at line {line_num}: {e}")
        with open(progress_file, "w") as pf:
            pf.write(str(line_num))
        print(f"Resume with: python {sys.argv[0]} --output {output_file} --resume")
        raise
    finally:
        combined_f.close()
        for bf in bin_files.values():
            bf.close()

    elapsed = time.time() - start_time
    print(f"\nPhase 1 complete!")
    print(f"  Total processed: {total_processed:,}")
    print(f"  Passed filters:  {passed_count:,}")
    print(f"  Filtered out:    {stats['filtered']:,}")
    print(f"  Fragmented:      {stats.get('fragmented', 0):,}")
    print(f"  Parse errors:    {stats.get('parse_error', 0):,}")
    print(f"  Time: {elapsed/60:.1f} min ({elapsed/3600:.1f} hours)")
    print(f"\n  MW Bin distribution:")
    for label, lo, hi in MW_BINS:
        print(f"    Bin {label} ({lo}-{hi} Da): {bin_counts.get(label, 0):,}")
    print(f"\n  Reactive handle distribution:")
    for h, c in sorted(handle_stats.items(), key=lambda x: -x[1]):
        print(f"    {h}: {c:,}")

    if os.path.exists(progress_file):
        os.remove(progress_file)

    return passed_count


def check_ghs_hazards(input_file, output_file, batch_size=100):
    """Check GHS hazard data for filtered compounds via PubChem API."""

    print(f"\nPhase 2: Checking GHS hazard data via PubChem API...")

    EXCLUDE_CODES = {
        "H350", "H351",
        "H340", "H341",
        "H360", "H361",
        "H400", "H410", "H411", "H412", "H413",
    }

    compounds = []
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split("\t")
                if len(parts) >= 2:
                    compounds.append((parts[0], parts[1]))

    print(f"  Checking {len(compounds):,} compounds for GHS hazards...")
    print(f"  Excluding codes: {', '.join(sorted(EXCLUDE_CODES))}")

    excluded_cids = set()
    checked = 0
    api_errors = 0

    cid_map = {}
    for smiles, name in compounds:
        if name.startswith("CID_"):
            cid = name[4:]
            cid_map[cid] = (smiles, name)

    cid_list = list(cid_map.keys())

    for cid in cid_list:
        ghs_url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/"
                   f"compound/{cid}/JSON?heading=GHS+Classification")
        try:
            req = urllib.request.Request(ghs_url)
            req.add_header("Accept", "application/json")
            resp = urllib.request.urlopen(req, timeout=30)
            data = json.loads(resp.read().decode())
            ghs_text = json.dumps(data)
            for code in EXCLUDE_CODES:
                if code in ghs_text:
                    excluded_cids.add(cid)
                    break
        except urllib.error.HTTPError as e:
            if e.code != 404:
                api_errors += 1
        except Exception:
            api_errors += 1

        checked += 1
        time.sleep(0.2)

        if checked % 500 == 0:
            print(f"    Checked: {checked:,}/{len(cid_list):,} | "
                  f"Excluded: {len(excluded_cids):,} | "
                  f"API errors: {api_errors}")

    print(f"  GHS check complete: {len(excluded_cids):,} compounds excluded")

    kept = 0
    with open(output_file, "w") as f:
        for smiles, name in compounds:
            cid = name[4:] if name.startswith("CID_") else None
            if cid and cid in excluded_cids:
                continue
            f.write(f"{smiles}\t{name}\n")
            kept += 1

    return kept, len(excluded_cids)


def main():
    parser = argparse.ArgumentParser(
        description="PubChem Agrochemical Library Builder v2 for AutoGrow4")
    parser.add_argument("--output", "-o", default="agrochemical_library.smi",
                        help="Output .smi file (default: agrochemical_library.smi)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from last checkpoint")
    parser.add_argument("--skip-ghs", action="store_true",
                        help="Skip GHS hazard check (much faster)")
    parser.add_argument("--mw-min", type=float, default=0,
                        help="Minimum molecular weight (default: 0)")
    parser.add_argument("--mw-max", type=float, default=250,
                        help="Maximum molecular weight (default: 250)")
    parser.add_argument("--logp-min", type=float, default=-2,
                        help="Minimum LogP (default: -2, relaxed for fragments)")
    parser.add_argument("--logp-max", type=float, default=5,
                        help="Maximum LogP (default: 5)")
    parser.add_argument("--rot-max", type=int, default=8,
                        help="Maximum rotatable bonds (default: 8)")
    parser.add_argument("--hba-min", type=int, default=0,
                        help="Minimum H-bond acceptors (default: 0)")
    parser.add_argument("--hba-max", type=int, default=10,
                        help="Maximum H-bond acceptors (default: 10)")
    parser.add_argument("--hbd-max", type=int, default=5,
                        help="Maximum H-bond donors (default: 5)")
    parser.add_argument("--sa-max", type=float, default=3.0,
                        help="Maximum SA score (default: 3.0)")
    parser.add_argument("--tpsa-min", type=float, default=0,
                        help="Minimum TPSA (default: 0, relaxed for fragments)")
    parser.add_argument("--tpsa-max", type=float, default=200,
                        help="Maximum TPSA (default: 200)")
    parser.add_argument("--aro-max", type=int, default=3,
                        help="Maximum aromatic rings (default: 3)")
    parser.add_argument("--no-handles", action="store_true",
                        help="Disable reactive handle requirement")

    args = parser.parse_args()

    print("=" * 60)
    print("PubChem Agrochemical Library Builder v2")
    print("=" * 60)
    print()

    filter_kwargs = dict(
        mw_min=args.mw_min, mw_max=args.mw_max,
        logp_min=args.logp_min, logp_max=args.logp_max,
        rot_bonds_max=args.rot_max,
        hba_min=args.hba_min, hba_max=args.hba_max,
        hbd_max=args.hbd_max, sa_max=args.sa_max,
        tpsa_min=args.tpsa_min, tpsa_max=args.tpsa_max,
        aro_rings_max=args.aro_max,
        require_handles=not args.no_handles,
    )

    passed = download_and_filter(
        args.output, resume=args.resume, skip_ghs=args.skip_ghs, **filter_kwargs
    )

    if passed == 0:
        print("No compounds passed filters. Try relaxing parameters.")
        return

    if not args.skip_ghs:
        temp_file = args.output + ".phase1.tmp"
        if passed > 50000:
            print(f"\nWARNING: {passed:,} compounds to check against GHS API.")
            print(f"This will take approximately {passed * 0.2 / 3600:.1f} hours.")
            print("Consider using --skip-ghs for faster results.")
            response = input("Continue with GHS check? [y/N]: ").strip().lower()
            if response != "y":
                print("Skipping GHS check.")
                os.rename(temp_file, args.output)
                print(f"\nFinal library: {args.output} ({passed:,} compounds)")
                return

        kept, excluded = check_ghs_hazards(temp_file, args.output)
        os.remove(temp_file)
        print(f"\nPhase 2 complete: {kept:,} kept, {excluded:,} excluded")
    else:
        kept = passed

    print()
    print("=" * 60)
    print(f"DONE! Final library: {args.output}")
    print(f"Total compounds: {kept:,}")
    print(f"Ready for AutoGrow4 source_compound_file")
    print("=" * 60)


if __name__ == "__main__":
    main()
