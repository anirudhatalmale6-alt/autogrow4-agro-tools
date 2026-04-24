#!/usr/bin/env python3
"""
PubChem Agrochemical Library Builder for AutoGrow4
===================================================
Downloads PubChem compound SMILES and filters for agrochemical drug-likeness.

Phase 1: Download CID-SMILES.gz from PubChem FTP, stream-filter with RDKit
Phase 2: Check GHS hazard data for filtered compounds via PubChem API
Phase 3: Output clean .smi file for AutoGrow4

Filters (Tice's Rules + custom):
  - MW 200-350
  - Elements: C, H, O, N, S only
  - LogP 0 to 4
  - Rotatable bonds < 12
  - H-bond acceptors 2-10
  - H-bond donors < 3
  - SA Score <= 3.0

Usage:
  conda activate ag4_full
  python pubchem_agro_library.py --output agrochemical_library.smi

  # Resume from where you left off (skips already-processed lines):
  python pubchem_agro_library.py --output agrochemical_library.smi --resume

  # Skip GHS hazard check:
  python pubchem_agro_library.py --output agrochemical_library.smi --skip-ghs
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
    from rdkit.Chem import Descriptors, Lipinski, Crippen
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    print("ERROR: RDKit not found. Activate conda environment first:")
    print("  conda activate ag4_full")
    sys.exit(1)

# SA Score setup
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
        print("WARNING: SA_Score module not found in RDKit Contrib.")
        print("Will use a simplified synthetic accessibility estimate.")
        print("For best results, install sascorer.py from RDKit Contrib/SA_Score/")
    except Exception:
        pass


PUBCHEM_FTP_URL = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz"
ALLOWED_ELEMENTS = {"C", "H", "O", "N", "S"}
CHUNK_REPORT = 500000


def calculate_sa_score(mol):
    """Calculate synthetic accessibility score (1=easy, 10=hard)."""
    if HAS_SASCORER:
        return sascorer.calculateScore(mol)
    # Fallback: rough estimate based on complexity metrics
    ring_count = Descriptors.RingCount(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    stereo_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    # Simple heuristic: more rings, stereo centers, and heavy atoms = harder
    score = 1.0 + (ring_count * 0.5) + (stereo_centers * 0.8) + (heavy_atoms * 0.05)
    return min(score, 10.0)


def get_elements(mol):
    """Get set of element symbols in a molecule."""
    return {atom.GetSymbol() for atom in mol.GetAtoms()}


def passes_filters(smiles, mw_min=200, mw_max=350, logp_min=0, logp_max=4,
                   rot_bonds_max=12, hba_min=2, hba_max=10, hbd_max=3,
                   sa_max=3.0):
    """Apply all filters to a SMILES string. Returns (pass, mol, props) or (False, None, None)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, None, None

    # Element check (fastest, do first)
    elements = get_elements(mol)
    if not elements.issubset(ALLOWED_ELEMENTS):
        return False, None, None

    # MW check
    mw = Descriptors.ExactMolWt(mol)
    if mw < mw_min or mw > mw_max:
        return False, None, None

    # LogP check (Tice's rule)
    logp = Crippen.MolLogP(mol)
    if logp < logp_min or logp > logp_max:
        return False, None, None

    # H-bond donors (Tice: < 3)
    hbd = Lipinski.NumHDonors(mol)
    if hbd >= hbd_max:
        return False, None, None

    # H-bond acceptors (Tice: 2-10)
    hba = Lipinski.NumHAcceptors(mol)
    if hba < hba_min or hba > hba_max:
        return False, None, None

    # Rotatable bonds (Tice: < 12)
    rot = Lipinski.NumRotatableBonds(mol)
    if rot >= rot_bonds_max:
        return False, None, None

    # SA Score (1=easy, 10=hard, want <= 3)
    sa = calculate_sa_score(mol)
    if sa > sa_max:
        return False, None, None

    props = {
        "mw": round(mw, 2),
        "logp": round(logp, 2),
        "hbd": hbd,
        "hba": hba,
        "rot": rot,
        "sa": round(sa, 2)
    }
    return True, mol, props


def download_and_filter(output_file, resume=False, mw_min=200, mw_max=350,
                        logp_min=0, logp_max=4, rot_bonds_max=12,
                        hba_min=2, hba_max=10, hbd_max=3, sa_max=3.0):
    """Download PubChem CID-SMILES.gz and filter compounds."""

    temp_file = output_file + ".phase1.tmp"
    progress_file = output_file + ".progress"
    skip_lines = 0

    if resume and os.path.exists(progress_file):
        with open(progress_file) as f:
            skip_lines = int(f.read().strip())
        print(f"Resuming from line {skip_lines}")

    mode = "a" if resume and os.path.exists(temp_file) else "w"

    stats = Counter()
    passed_count = 0
    total_processed = 0

    print(f"Downloading PubChem CID-SMILES from: {PUBCHEM_FTP_URL}")
    print(f"Filters: MW {mw_min}-{mw_max}, LogP {logp_min}-{logp_max}, "
          f"HBD <{hbd_max}, HBA {hba_min}-{hba_max}, RotBonds <{rot_bonds_max}, "
          f"SA <={sa_max}, Elements: {','.join(sorted(ALLOWED_ELEMENTS))}")
    print(f"Output: {temp_file}")
    print(f"SA Score method: {'sascorer (accurate)' if HAS_SASCORER else 'heuristic (approximate)'}")
    print()

    start_time = time.time()

    try:
        req = urllib.request.Request(PUBCHEM_FTP_URL)
        response = urllib.request.urlopen(req, timeout=120)
        # Wrap in gzip decompressor
        gz = gzip.GzipFile(fileobj=response)

        with open(temp_file, mode) as out_f:
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

                passed, mol, props = passes_filters(
                    smiles, mw_min, mw_max, logp_min, logp_max,
                    rot_bonds_max, hba_min, hba_max, hbd_max, sa_max
                )

                if passed:
                    out_f.write(f"{smiles}\tCID_{cid}\n")
                    passed_count += 1
                    stats["passed"] += 1
                else:
                    stats["filtered"] += 1

                if total_processed % CHUNK_REPORT == 0:
                    elapsed = time.time() - start_time
                    rate = total_processed / elapsed if elapsed > 0 else 0
                    print(f"  Processed: {total_processed:>12,} | "
                          f"Passed: {passed_count:>8,} | "
                          f"Rate: {rate:,.0f}/sec | "
                          f"Elapsed: {elapsed/60:.1f} min")
                    # Save progress
                    with open(progress_file, "w") as pf:
                        pf.write(str(line_num))
                    out_f.flush()

        gz.close()
        response.close()

    except KeyboardInterrupt:
        print(f"\nInterrupted at line {line_num}. Progress saved.")
        with open(progress_file, "w") as pf:
            pf.write(str(line_num))
        print(f"Resume with: python {sys.argv[0]} --output {output_file} --resume")
        return temp_file, passed_count

    except Exception as e:
        print(f"\nError at line {line_num}: {e}")
        with open(progress_file, "w") as pf:
            pf.write(str(line_num))
        print(f"Resume with: python {sys.argv[0]} --output {output_file} --resume")
        raise

    elapsed = time.time() - start_time
    print(f"\nPhase 1 complete!")
    print(f"  Total processed: {total_processed:,}")
    print(f"  Passed filters:  {passed_count:,}")
    print(f"  Filtered out:    {stats['filtered']:,}")
    print(f"  Parse errors:    {stats.get('parse_error', 0):,}")
    print(f"  Time: {elapsed/60:.1f} min ({elapsed/3600:.1f} hours)")

    # Clean up progress file
    if os.path.exists(progress_file):
        os.remove(progress_file)

    return temp_file, passed_count


def check_ghs_hazards(input_file, output_file, batch_size=100):
    """Check GHS hazard data for filtered compounds via PubChem API.
    Removes compounds flagged as carcinogenic, mutagenic, or environmentally toxic."""

    print(f"\nPhase 2: Checking GHS hazard data via PubChem API...")

    # GHS hazard codes to exclude
    # H350/H351: Carcinogenic
    # H340/H341: Mutagenic
    # H360/H361: Reproductive toxicity
    # H400/H410/H411: Environmental toxicity (aquatic)
    EXCLUDE_CODES = {
        "H350", "H351",  # Carcinogenicity
        "H340", "H341",  # Mutagenicity
        "H360", "H361",  # Reproductive toxicity
        "H400", "H410", "H411", "H412", "H413",  # Aquatic toxicity
    }

    # Read all compounds from phase 1
    compounds = []
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split("\t")
                if len(parts) >= 2:
                    compounds.append((parts[0], parts[1]))  # smiles, name

    print(f"  Checking {len(compounds):,} compounds for GHS hazards...")
    print(f"  Excluding codes: {', '.join(sorted(EXCLUDE_CODES))}")

    excluded_cids = set()
    checked = 0
    api_errors = 0

    # Extract CIDs for API queries
    cid_map = {}
    for smiles, name in compounds:
        if name.startswith("CID_"):
            cid = name[4:]
            cid_map[cid] = (smiles, name)

    cid_list = list(cid_map.keys())

    # Query in batches
    for i in range(0, len(cid_list), batch_size):
        batch = cid_list[i:i + batch_size]
        cids_str = ",".join(batch)

        url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
               f"{cids_str}/property/CanonicalSMILES/JSON")

        # Use PUG-VIEW for GHS data
        for cid in batch:
            ghs_url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/"
                       f"compound/{cid}/JSON?heading=GHS+Classification")
            try:
                req = urllib.request.Request(ghs_url)
                req.add_header("Accept", "application/json")
                resp = urllib.request.urlopen(req, timeout=30)
                data = json.loads(resp.read().decode())

                # Parse GHS data for hazard codes
                ghs_text = json.dumps(data)
                for code in EXCLUDE_CODES:
                    if code in ghs_text:
                        excluded_cids.add(cid)
                        break

            except urllib.error.HTTPError as e:
                if e.code == 404:
                    pass  # No GHS data = keep the compound
                else:
                    api_errors += 1
            except Exception:
                api_errors += 1

            checked += 1
            time.sleep(0.2)  # Rate limit: 5 req/sec

            if checked % 500 == 0:
                print(f"    Checked: {checked:,}/{len(cid_list):,} | "
                      f"Excluded: {len(excluded_cids):,} | "
                      f"API errors: {api_errors}")

    print(f"  GHS check complete: {len(excluded_cids):,} compounds excluded")

    # Write final output
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
        description="PubChem Agrochemical Library Builder for AutoGrow4")
    parser.add_argument("--output", "-o", default="agrochemical_library.smi",
                        help="Output .smi file (default: agrochemical_library.smi)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from last checkpoint")
    parser.add_argument("--skip-ghs", action="store_true",
                        help="Skip GHS hazard check (much faster)")
    parser.add_argument("--mw-min", type=float, default=200,
                        help="Minimum molecular weight (default: 200)")
    parser.add_argument("--mw-max", type=float, default=350,
                        help="Maximum molecular weight (default: 350)")
    parser.add_argument("--logp-min", type=float, default=0,
                        help="Minimum LogP (default: 0)")
    parser.add_argument("--logp-max", type=float, default=4,
                        help="Maximum LogP (default: 4)")
    parser.add_argument("--rot-max", type=int, default=12,
                        help="Maximum rotatable bonds (default: 12)")
    parser.add_argument("--hba-min", type=int, default=2,
                        help="Minimum H-bond acceptors (default: 2)")
    parser.add_argument("--hba-max", type=int, default=10,
                        help="Maximum H-bond acceptors (default: 10)")
    parser.add_argument("--hbd-max", type=int, default=3,
                        help="Maximum H-bond donors (default: 3)")
    parser.add_argument("--sa-max", type=float, default=3.0,
                        help="Maximum SA score (default: 3.0)")

    args = parser.parse_args()

    print("=" * 60)
    print("PubChem Agrochemical Library Builder")
    print("=" * 60)
    print()

    # Phase 1: Download and filter
    temp_file, passed = download_and_filter(
        args.output, resume=args.resume,
        mw_min=args.mw_min, mw_max=args.mw_max,
        logp_min=args.logp_min, logp_max=args.logp_max,
        rot_bonds_max=args.rot_max,
        hba_min=args.hba_min, hba_max=args.hba_max,
        hbd_max=args.hbd_max, sa_max=args.sa_max
    )

    if passed == 0:
        print("No compounds passed filters. Try relaxing parameters.")
        return

    # Phase 2: GHS hazard check (optional)
    if not args.skip_ghs:
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
        os.rename(temp_file, args.output)
        kept = passed

    print()
    print("=" * 60)
    print(f"DONE! Final library: {args.output}")
    print(f"Total compounds: {kept:,}")
    print(f"Ready for AutoGrow4 source_compound_file")
    print("=" * 60)


if __name__ == "__main__":
    main()
