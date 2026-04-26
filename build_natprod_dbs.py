#!/usr/bin/env python3
"""
Download and build natural product databases on Expanse.
Creates np_curated.smi and np_huggingface.smi in ~/final_project/natural_products/
"""

import csv
import gzip
import io
import os
import sys
import urllib.request

OUT_DIR = os.path.expanduser("~/final_project/natural_products")
os.makedirs(OUT_DIR, exist_ok=True)

AGONISTS_DECOYS = [
    # Agonists - SMILES from PubChem verified entries
    ("CC12C(CCC3(C1C(C45C3CCC(C4)C(=C)C5)C(=O)O)OC2=O)O", "GA4_agonist"),  # CID 92109
    ("CC12C(C=CC3(C1C(C45C3CCC(C4)(C(=C)C5)O)C(=O)O)OC2=O)O", "GA3_agonist"),  # CID 6466
    ("CC12C(C=CC3(C1C(C45C3CCC(C4)C(=C)C5)C(=O)O)OC2=O)O", "GA7_agonist"),  # CID 92782
    ("CC12C=CCC3(C1C(C45C3CCC(C4)(C(=C)C5)O)C(=O)O)OC2=O", "GA5_agonist"),  # CID 5281988
    ("C1CCC(CC1)(C(=O)N)N2C(=O)C3=C(C2=O)C(=CC=C3)Cl", "AC94377_agonist"),  # CID 93267
    ("CC1=CC=CC=C1NS(=O)(=O)C2=CC=C(C=C2)C(=O)O", "Chemical6_agonist"),
    ("O=C1C2=CC=CC=C2C(=O)N1CC3=CC=CC=C3", "Chemical2_agonist"),
    ("CC1=CC=C(C=C1)N2C(=O)C3=CC=CC=C3C2=O", "Chemical5_agonist"),
    ("O=C(NC1=CC=C(C=C1)C(=O)O)C2=CC=C(C=C2)Cl", "Chemical8_agonist"),
    # Decoys - SMILES from PubChem/KEGG verified entries
    ("CC12C3C(C45CC(=C)C(C4)(CCC5C3(CC(C1O)O)OC2=O)O)C(=O)O", "GA8_decoy"),  # CID 5280607
    ("CC12CC(CC3(C1C(C45C3CCC(C4)(C(=C)C5)O)C(=O)O)OC2=O)O", "GA29_decoy"),  # CID 14605548
    ("CC12C3C(C45CC(CCC4C3(CC(C1O)O)OC2=O)C(=C)C5)C(=O)O", "GA34_decoy"),  # CID 5281987
    ("CC12CCCC(C1C(C34C2CCC(C3)C(=C)C4)C(=O)O)(C)C(=O)O", "GA12_decoy"),  # CID 443450
    ("C=C1C[C@]23C[C@H]1CC[C@H]2[C@@]12CCC[C@@](C)(C(=O)OC1)[C@H]2[C@@H]3C(=O)O", "GA15_decoy"),  # KEGG C14162
    ("CC12CCCC(C1C(C34C2CCC(C3)(C(=C)C4)O)C(=O)O)(C)C(=O)O", "GA53_decoy"),  # CID 440914
    ("O=C1NC(=O)C2=CC=CC=C12", "AC94377Analog_decoy"),
]


def download(url, desc=""):
    print(f"  Downloading {desc or url}...")
    sys.stdout.flush()
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    resp = urllib.request.urlopen(req, timeout=600)
    data = resp.read()
    print(f"  Downloaded {len(data)/1024/1024:.1f} MB")
    sys.stdout.flush()
    return data


def build_npatlas():
    print("\n[1/5] Natural Products Atlas...")
    url = "https://www.npatlas.org/static/downloads/NPAtlas_download.tsv"
    data = download(url, "NP Atlas TSV")
    reader = csv.DictReader(io.StringIO(data.decode("utf-8")), delimiter="\t")
    compounds = []
    for row in reader:
        smi = row.get("compound_smiles", "").strip()
        name = row.get("compound_names", row.get("npaid", "")).strip()
        if smi and name:
            name = name.replace("\t", "_").replace(" ", "_")
            compounds.append((smi, name))
    print(f"  NP Atlas: {len(compounds)} compounds")
    return compounds


def build_chebi():
    print("\n[2/5] ChEBI...")
    try:
        url = "https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/structures.csv.gz"
        data = download(url, "ChEBI structures")
        struct_data = gzip.decompress(data).decode("utf-8", errors="replace")

        url_names = "https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv.gz"
        data2 = download(url_names, "ChEBI compounds")
        names_data = gzip.decompress(data2).decode("utf-8", errors="replace")

        id_to_name = {}
        for line in names_data.split("\n")[1:]:
            parts = line.split("\t")
            if len(parts) >= 6:
                cid = parts[0].strip()
                name = parts[5].strip() if parts[5].strip() else parts[0].strip()
                name = name.replace("\t", "_").replace(" ", "_")
                id_to_name[cid] = name

        compounds = []
        seen = set()
        reader = csv.reader(io.StringIO(struct_data))
        header = next(reader, None)
        if header:
            for row in reader:
                if len(row) >= 4:
                    cid = row[0].strip()
                    struct_type = row[1].strip()
                    structure = row[2].strip()
                    if struct_type == "SMILES" and structure:
                        if structure not in seen:
                            seen.add(structure)
                            name = id_to_name.get(cid, f"CHEBI_{cid}")
                            compounds.append((structure, name))

        print(f"  ChEBI: {len(compounds)} compounds")
        return compounds
    except Exception as e:
        print(f"  ChEBI failed: {e}")
        print("  Trying alternative ChEBI download...")
        try:
            url2 = "https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chebiId_inchi.tsv"
            data = download(url2, "ChEBI InChI TSV")
            # Fallback: use ChEBI SDF lite
            url3 = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete_3star.sdf.gz"
            data = download(url3, "ChEBI SDF (3-star)")
            sdf_text = gzip.decompress(data).decode("utf-8", errors="replace")
            compounds = []
            seen = set()
            from rdkit import Chem
            suppl = Chem.SDMolSupplier()
            for block in sdf_text.split("$$$$\n"):
                if not block.strip():
                    continue
                try:
                    mol = Chem.MolFromMolBlock(block)
                    if mol:
                        smi = Chem.MolToSmiles(mol)
                        name = mol.GetProp("ChEBI Name") if mol.HasProp("ChEBI Name") else f"CHEBI_{len(compounds)}"
                        name = name.replace("\t", "_").replace(" ", "_")
                        if smi not in seen:
                            seen.add(smi)
                            compounds.append((smi, name))
                except Exception:
                    pass
            print(f"  ChEBI (SDF fallback): {len(compounds)} compounds")
            return compounds
        except Exception as e2:
            print(f"  ChEBI fallback also failed: {e2}")
            return []


def build_phytochemdb():
    print("\n[3/5] PhytochemDB...")
    try:
        url = "https://www.phytochemdb.com/static/SDF_Phytochemdb.sdf"
        data = download(url, "PhytochemDB SDF")
        sdf_text = data.decode("utf-8", errors="replace")

        compounds = []
        blocks = sdf_text.split("$$$$")
        for i, block in enumerate(blocks):
            block = block.strip()
            if not block:
                continue
            lines = block.split("\n")
            name = lines[0].strip() if lines else f"Phyto_{i}"
            name = name.replace("\t", "_").replace(" ", "_") or f"Phyto_{i}"

            # Try to find SMILES in properties
            smi = None
            for j, line in enumerate(lines):
                if "SMILES" in line.upper() and j + 1 < len(lines):
                    candidate = lines[j + 1].strip()
                    if candidate and not candidate.startswith(">"):
                        smi = candidate
                        break

            if not smi:
                try:
                    from rdkit import Chem
                    mol = Chem.MolFromMolBlock(block + "\n$$$$")
                    if mol:
                        smi = Chem.MolToSmiles(mol)
                except Exception:
                    pass

            if smi:
                compounds.append((smi, name))

        print(f"  PhytochemDB: {len(compounds)} compounds")
        return compounds
    except Exception as e:
        print(f"  PhytochemDB failed: {e}")
        return []


def build_biogem():
    print("\n[4/5] BioGEM PDTDB...")
    compounds = []
    for i in range(1, 238):
        lid = f"pdtdbl{i:05d}"
        url = f"https://pdt.biogem.org/output.php?ligand_id={lid}"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            resp = urllib.request.urlopen(req, timeout=30)
            html = resp.read().decode("utf-8", errors="replace")

            smi = None
            name = lid
            for line in html.split("\n"):
                if "Canonical SMILES" in line:
                    import re
                    match = re.search(r'<td[^>]*>([^<]+)</td>', html[html.index("Canonical SMILES"):])
                    if match:
                        candidate = match.group(1).strip()
                        if candidate and len(candidate) > 2:
                            smi = candidate

                if "Ligand Name" in line:
                    match = re.search(r'Ligand Name.*?<td[^>]*>([^<]+)</td>', html, re.DOTALL)
                    if match:
                        name = match.group(1).strip().replace("\t", "_").replace(" ", "_")

            if smi:
                compounds.append((smi, name))
        except Exception:
            pass

        if i % 50 == 0:
            print(f"  BioGEM: {i}/237 processed, {len(compounds)} found...")
            sys.stdout.flush()

    print(f"  BioGEM PDTDB: {len(compounds)} compounds")
    return compounds


def build_huggingface():
    print("\n[5/5] HuggingFace Bioactives/Naturals...")
    url = "https://huggingface.co/datasets/gbyuvd/bioactives-naturals-smiles-molgen/resolve/main/comb_smi.csv"
    data = download(url, "HF Bioactives CSV (~154MB)")
    lines = data.decode("utf-8", errors="replace").split("\n")

    compounds = []
    seen = set()
    for i, line in enumerate(lines):
        smi = line.strip()
        if not smi or smi == "smiles" or smi == "SMILES":
            continue
        if smi not in seen:
            seen.add(smi)
            compounds.append((smi, f"HF_BIOACTIVE_{i}"))

        if (i + 1) % 500000 == 0:
            print(f"  Processed {i+1:,} lines, {len(compounds):,} unique...")
            sys.stdout.flush()

    print(f"  HuggingFace: {len(compounds)} unique compounds")
    return compounds


def main():
    print("=" * 60)
    print("Natural Products Database Builder for Vina Screening")
    print("=" * 60)

    # Build curated databases (each wrapped in try/except so one failure doesn't kill all)
    npatlas = []
    chebi = []
    phytochem = []
    biogem = []

    try:
        npatlas = build_npatlas()
    except Exception as e:
        print(f"  NP Atlas FAILED: {e}")

    try:
        chebi = build_chebi()
    except Exception as e:
        print(f"  ChEBI FAILED: {e}")

    try:
        phytochem = build_phytochemdb()
    except Exception as e:
        print(f"  PhytochemDB FAILED: {e}")

    try:
        biogem = build_biogem()
    except Exception as e:
        print(f"  BioGEM FAILED: {e}")

    # Write curated file
    curated_path = os.path.join(OUT_DIR, "np_curated.smi")
    seen = set()
    count = 0
    with open(curated_path, "w") as f:
        for smi, name in AGONISTS_DECOYS:
            f.write(f"{smi}\tCONTROLS::{name}\n")
            seen.add(smi)
            count += 1

        for compounds, prefix in [(biogem, "BioGEM_PDTDB"),
                                   (npatlas, "NP_Atlas"),
                                   (phytochem, "PhytochemDB"),
                                   (chebi, "ChEBI")]:
            added = 0
            for smi, name in compounds:
                if smi not in seen:
                    seen.add(smi)
                    f.write(f"{smi}\t{prefix}::{name}\n")
                    count += 1
                    added += 1
            print(f"  {prefix}: {added} unique added to curated")

    print(f"\nCurated total: {count} compounds -> {curated_path}")

    # Build HuggingFace
    hf_compounds = []
    try:
        hf_compounds = build_huggingface()
    except Exception as e:
        print(f"  HuggingFace FAILED: {e}")

    hf_path = os.path.join(OUT_DIR, "np_huggingface.smi")
    seen_hf = set()
    count_hf = 0
    with open(hf_path, "w") as f:
        for smi, name in AGONISTS_DECOYS:
            f.write(f"{smi}\tCONTROLS::{name}\n")
            seen_hf.add(smi)
            count_hf += 1

        for smi, name in hf_compounds:
            if smi not in seen_hf:
                seen_hf.add(smi)
                f.write(f"{smi}\tHF_Bioactives::{name}\n")
                count_hf += 1

    print(f"HuggingFace total: {count_hf} compounds -> {hf_path}")
    print("\nDone! Files ready for screening.")


if __name__ == "__main__":
    main()
