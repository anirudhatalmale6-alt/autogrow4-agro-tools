#!/usr/bin/env python3
"""
Water-Validated Cross-Docking - Task 6

Re-runs the cross-docking IFP matrix from Task 3 with conserved
water molecules embedded in the receptor. Compares "dry" vs "wet"
docking to validate that the placed waters improve:
  1. CNN scores (more favorable interactions)
  2. IFP Tanimoto (water-mediated contacts captured)
  3. RMSD (waters constrain ligand orientation)

Water sets tested:
  - Dry: standard receptor (Task 3 baseline)
  - Wet: receptor + pharmacophore waters + Asp289 shell
"""

import os
import sys
import json
import subprocess
import urllib.request
import numpy as np

from Bio.PDB import PDBParser, Superimposer, PDBIO, Select

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

OUTPUT_DIR = os.path.expanduser("~/final_project/water_validated_docking")
POCKET_DIR = os.path.expanduser("~/final_project/pocket_analysis")
TASK3_DIR = os.path.expanduser("~/final_project/cross_docking")
WATER_DIR = os.path.expanduser("~/final_project/water_mapping")
GNINA_PATH = os.path.expanduser("~/final_project/gnina")
os.makedirs(OUTPUT_DIR, exist_ok=True)

PDB_IDS = ["2ZSH", "2ZSI", "3ED1", "3EBL"]
SPECIES = {"2ZSH": "Arabidopsis", "2ZSI": "Arabidopsis",
           "3ED1": "Rice", "3EBL": "Rice"}
LIGAND_OF = {"2ZSH": "GA3", "2ZSI": "GA4",
             "3ED1": "GA3", "3EBL": "GA4"}

ARAB_RESIDUES = [24, 27, 31, 116, 119, 126, 127, 191, 218, 238, 239, 244, 319, 323]
RICE_RESIDUES = [24, 27, 31, 123, 126, 133, 134, 198, 225, 245, 246, 251, 326, 330]
BACKBONE_ATOMS = ["N", "CA", "C", "O"]

ARAB_TO_RICE = {
    24: 24, 27: 27, 31: 31, 116: 123, 119: 126, 126: 133, 127: 134,
    191: 198, 218: 225, 238: 245, 239: 246, 244: 251, 289: 296,
    319: 326, 323: 330,
}

IFP_TYPES = ["hbond_don", "hbond_acc", "hydrophobic", "aromatic",
             "ionic", "water_mediated"]
IFP_CUTOFFS = {
    "hbond_don": 3.5, "hbond_acc": 3.5, "hydrophobic": 4.5,
    "aromatic": 5.5, "ionic": 4.0, "water_mediated": 4.0,
}

HBOND_DONORS = {"ARG", "LYS", "HIS", "ASN", "GLN", "SER", "THR",
                "TRP", "TYR", "CYS"}
HBOND_ACCEPTORS = {"ASP", "GLU", "ASN", "GLN", "SER", "THR", "HIS",
                   "TYR", "CYS"}
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "TRP", "MET"}
AROMATIC_RES = {"PHE", "TYR", "TRP", "HIS"}
IONIC_POS = {"ARG", "LYS", "HIS"}
IONIC_NEG = {"ASP", "GLU"}

EXHAUSTIVENESS = 16
TC_THRESHOLD = 0.85


def download_pdb(pdb_id):
    path = os.path.join(POCKET_DIR, f"{pdb_id}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        return path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    urllib.request.urlretrieve(url, path)
    return path


def get_pocket_residues(pdb_id):
    return list(ARAB_RESIDUES) if SPECIES[pdb_id] == "Arabidopsis" \
        else list(RICE_RESIDUES)


def map_residue(arab_resnum, pdb_id):
    if SPECIES[pdb_id] == "Arabidopsis":
        return arab_resnum
    return ARAB_TO_RICE.get(arab_resnum, arab_resnum)


def get_pocket_backbone_atoms(chain, pocket_residues):
    atoms = []
    for resnum in pocket_residues:
        try:
            res = chain[(' ', resnum, ' ')]
            for aname in BACKBONE_ATOMS:
                if aname in res:
                    atoms.append(res[aname])
        except KeyError:
            pass
    return atoms


class ChainAProteinSelect(Select):
    def accept_chain(self, chain):
        return chain.id == "A"
    def accept_residue(self, residue):
        return residue.id[0] == " "


def extract_receptor(struct, pdb_id):
    out_path = os.path.join(OUTPUT_DIR, f"{pdb_id}_receptor_dry.pdb")
    io = PDBIO()
    io.set_structure(struct)
    io.save(out_path, ChainAProteinSelect())
    return out_path


def extract_aligned_ligand(struct, pdb_id, lig_name):
    out_path = os.path.join(OUTPUT_DIR,
                            f"{pdb_id}_{lig_name}_ligand.pdb")
    chain = struct[0]['A']
    lines = []
    serial = 1
    for res in chain:
        if res.get_resname() == lig_name and res.id[0] != " ":
            for atom in res:
                coord = atom.get_vector().get_array()
                elem = atom.element.strip()
                if not elem:
                    elem = atom.name.strip()[0]
                line = (f"HETATM{serial:5d} {atom.name:<4s} "
                        f"{lig_name:>3s} A{res.id[1]:4d}    "
                        f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                        f"  1.00  0.00          {elem:>2s}\n")
                lines.append(line)
                serial += 1
    if not lines:
        return None
    with open(out_path, "w") as f:
        for line in lines:
            f.write(line)
        f.write("END\n")
    return out_path


def get_pharmacophore_waters(struct, pdb_id):
    chain = struct[0]['A']
    lig_name = LIGAND_OF[pdb_id]
    waters_to_embed = []

    pharm_atoms = {
        "7-COOH": ["O71", "O72"],
        "3-OH": ["O31"],
        "13-OH": ["O13"],
    }

    asp_resnum = map_residue(289, pdb_id)

    lig_res = None
    for res in chain:
        if res.get_resname() == lig_name and res.id[0] != " ":
            lig_res = res
            break

    all_waters = []
    for res in chain:
        if res.get_resname() == "HOH":
            for a in res.get_atoms():
                if a.element == "O":
                    all_waters.append((res, a))

    seen_ids = set()

    for pharm_name, atom_names in pharm_atoms.items():
        if pharm_name == "13-OH" and lig_name == "GA4":
            continue
        if lig_res is None:
            continue

        for aname in atom_names:
            if aname not in lig_res:
                continue
            pharm_atom = lig_res[aname]
            pharm_coord = pharm_atom.get_vector().get_array()

            for wat_res, wat_o in all_waters:
                if wat_res.id[1] in seen_ids:
                    continue
                dist = np.linalg.norm(
                    wat_o.get_vector().get_array() - pharm_coord)
                if dist <= 3.5:
                    waters_to_embed.append({
                        "resid": wat_res.id[1],
                        "coord": wat_o.get_vector().get_array(),
                        "bfactor": wat_o.get_bfactor(),
                        "source": pharm_name,
                    })
                    seen_ids.add(wat_res.id[1])

    try:
        asp_res = chain[(' ', asp_resnum, ' ')]
        asp_atoms = list(asp_res.get_atoms())
        for wat_res, wat_o in all_waters:
            if wat_res.id[1] in seen_ids:
                continue
            min_d = min(
                np.linalg.norm(wat_o.get_vector().get_array() -
                               a.get_vector().get_array())
                for a in asp_atoms
            )
            if min_d <= 3.5:
                waters_to_embed.append({
                    "resid": wat_res.id[1],
                    "coord": wat_o.get_vector().get_array(),
                    "bfactor": wat_o.get_bfactor(),
                    "source": "Asp289_shell",
                })
                seen_ids.add(wat_res.id[1])
    except KeyError:
        pass

    return waters_to_embed


def create_wet_receptor(dry_receptor_path, waters, pdb_id):
    wet_path = os.path.join(OUTPUT_DIR, f"{pdb_id}_receptor_wet.pdb")

    with open(dry_receptor_path) as f:
        lines = f.readlines()

    max_serial = 0
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            try:
                s = int(line[6:11])
                if s > max_serial:
                    max_serial = s
            except ValueError:
                pass

    pdb_lines = [l for l in lines if not l.startswith(("END", "TER"))]

    serial = max_serial + 1
    for w in waters:
        coord = w["coord"]
        bf = w["bfactor"]
        resid = w["resid"]
        line = (f"HETATM{serial:5d}  O   HOH A{resid:4d}    "
                f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                f"  1.00{bf:6.2f}           O\n")
        pdb_lines.append(line)
        serial += 1

    pdb_lines.append("END\n")

    with open(wet_path, "w") as f:
        f.writelines(pdb_lines)

    return wet_path


def extract_ligand_sdf(lig_pdb_path, pdb_id, lig_name):
    out_path = os.path.join(OUTPUT_DIR, f"{pdb_id}_{lig_name}_ligand.sdf")
    if not HAS_RDKIT:
        return None
    mol = Chem.MolFromPDBFile(lig_pdb_path, removeHs=False, sanitize=False)
    if mol is None:
        mol = Chem.MolFromPDBFile(lig_pdb_path, removeHs=True, sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    writer = Chem.SDWriter(out_path)
    writer.write(mol)
    writer.close()
    return out_path


def run_gnina_dock(receptor, ligand_sdf, autobox_lig, out_sdf, label=""):
    if os.path.exists(out_sdf) and os.path.getsize(out_sdf) > 100:
        return True

    cmd = [
        GNINA_PATH,
        "-r", receptor,
        "-l", ligand_sdf,
        "-o", out_sdf,
        "--autobox_ligand", autobox_lig,
        "--autobox_add", "4",
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", "9",
        "--cnn_scoring", "rescore",
    ]

    cuda_devs = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if not cuda_devs:
        cmd.append("--no_gpu")

    try:
        env = os.environ.copy()
        env.setdefault("OMP_NUM_THREADS", "4")
        result = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=1800, env=env)
        if result.returncode != 0:
            print(f"      GNINA error: {result.stderr[:200]}")
            return False
        return os.path.exists(out_sdf) and os.path.getsize(out_sdf) > 50
    except Exception as e:
        print(f"      Exception: {e}")
        return False


def get_gnina_scores(docked_sdf):
    if not HAS_RDKIT or not os.path.exists(docked_sdf):
        return None, None
    suppl = Chem.SDMolSupplier(docked_sdf, removeHs=False)
    for mol in suppl:
        if mol is None:
            continue
        props = mol.GetPropsAsDict()
        return props.get("CNNscore"), props.get("minimizedAffinity")
    return None, None


def calc_rmsd(docked_sdf, crystal_pdb):
    if not HAS_RDKIT:
        return None
    ref = Chem.MolFromPDBFile(crystal_pdb, removeHs=True, sanitize=False)
    if ref is None:
        return None
    try:
        Chem.SanitizeMol(ref)
    except Exception:
        pass

    suppl = Chem.SDMolSupplier(docked_sdf, removeHs=True)
    best = 999.0
    for mol in suppl:
        if mol is None:
            continue
        try:
            n = min(ref.GetNumAtoms(), mol.GetNumAtoms())
            c1 = np.array([ref.GetConformer().GetAtomPosition(i)
                           for i in range(n)])
            c2 = np.array([mol.GetConformer().GetAtomPosition(i)
                           for i in range(n)])
            rmsd = np.sqrt(np.mean(np.sum((c1 - c2) ** 2, axis=1)))
            if rmsd < best:
                best = rmsd
        except Exception:
            continue
    return best if best < 900 else None


def compute_ifp(receptor_pdb, docked_sdf, pdb_id, include_waters=False):
    pocket_res = get_pocket_residues(pdb_id)
    n_types = len(IFP_TYPES)
    n_res = len(pocket_res)
    fp = np.zeros(n_types * n_res, dtype=int)

    res_atoms = {}
    water_coords = []
    with open(receptor_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    resnum = int(line[22:26])
                    if resnum in pocket_res:
                        resname = line[17:20].strip()
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        if resnum not in res_atoms:
                            res_atoms[resnum] = {"name": resname, "coords": []}
                        res_atoms[resnum]["coords"].append(np.array([x, y, z]))
                except (ValueError, IndexError):
                    pass
            elif include_waters and line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname == "HOH":
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        water_coords.append(np.array([x, y, z]))
                    except (ValueError, IndexError):
                        pass

    if not HAS_RDKIT or not os.path.exists(docked_sdf):
        return fp

    suppl = Chem.SDMolSupplier(docked_sdf, removeHs=True)
    lig_coords = []
    for mol in suppl:
        if mol is None:
            continue
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            lig_coords.append(np.array([pos.x, pos.y, pos.z]))
        break

    if not lig_coords:
        return fp

    lig_coords = np.array(lig_coords)

    for res_idx, resnum in enumerate(pocket_res):
        if resnum not in res_atoms:
            continue
        resname = res_atoms[resnum]["name"]
        rcoords = np.array(res_atoms[resnum]["coords"])

        dists = np.min([np.linalg.norm(rcoords - lc, axis=1)
                        for lc in lig_coords], axis=0)
        min_dist = np.min(dists)

        for type_idx, itype in enumerate(IFP_TYPES):
            bit_idx = type_idx * n_res + res_idx
            cutoff = IFP_CUTOFFS[itype]

            if itype == "water_mediated" and include_waters:
                for wc in water_coords:
                    d_res = np.min(np.linalg.norm(rcoords - wc, axis=1))
                    d_lig = np.min(np.linalg.norm(lig_coords - wc, axis=1))
                    if d_res <= cutoff and d_lig <= cutoff:
                        fp[bit_idx] = 1
                        break
            elif itype == "hbond_don" and min_dist <= cutoff:
                if resname in HBOND_DONORS:
                    fp[bit_idx] = 1
            elif itype == "hbond_acc" and min_dist <= cutoff:
                if resname in HBOND_ACCEPTORS:
                    fp[bit_idx] = 1
            elif itype == "hydrophobic" and min_dist <= cutoff:
                if resname in HYDROPHOBIC_RES:
                    fp[bit_idx] = 1
            elif itype == "aromatic" and min_dist <= cutoff:
                if resname in AROMATIC_RES:
                    fp[bit_idx] = 1
            elif itype == "ionic" and min_dist <= cutoff:
                if resname in (IONIC_POS | IONIC_NEG):
                    fp[bit_idx] = 1

    return fp


def tanimoto(fp1, fp2):
    a_and_b = np.sum(fp1 & fp2)
    a_or_b = np.sum(fp1 | fp2)
    if a_or_b == 0:
        return 0.0
    return float(a_and_b) / float(a_or_b)


def main():
    print("=" * 70)
    print("WATER-VALIDATED CROSS-DOCKING - Task 6")
    print("=" * 70)
    print(f"Output: {OUTPUT_DIR}")
    print(f"GNINA: {GNINA_PATH}")
    print(f"Exhaustiveness: {EXHAUSTIVENESS}")
    print(f"IFP types: {IFP_TYPES}")
    print(f"Tanimoto threshold: {TC_THRESHOLD}")
    print()

    has_gnina = os.path.isfile(GNINA_PATH)
    if not has_gnina:
        print("ERROR: GNINA not found. Cannot proceed.")
        print(f"  Expected: {GNINA_PATH}")
        sys.exit(1)

    parser = PDBParser(QUIET=True)
    structures = {}
    chains = {}

    print("Step 1: Loading and aligning structures...")
    for pdb_id in PDB_IDS:
        pdb_path = download_pdb(pdb_id)
        struct = parser.get_structure(pdb_id, pdb_path)
        structures[pdb_id] = struct
        chains[pdb_id] = struct[0]['A']

    ref_id = "2ZSH"
    ref_bb = get_pocket_backbone_atoms(
        chains[ref_id], get_pocket_residues(ref_id))

    for pdb_id in PDB_IDS:
        if pdb_id == ref_id:
            continue
        mov_bb = get_pocket_backbone_atoms(
            chains[pdb_id], get_pocket_residues(pdb_id))
        n = min(len(ref_bb), len(mov_bb))
        sup = Superimposer()
        sup.set_atoms(ref_bb[:n], mov_bb[:n])
        sup.apply(structures[pdb_id][0].get_atoms())
        print(f"  {pdb_id} aligned: RMSD = {sup.rms:.4f} A")

    print(f"\nStep 2: Preparing dry and wet receptors...")
    dry_receptors = {}
    wet_receptors = {}
    lig_pdbs = {}
    lig_sdfs = {}
    water_counts = {}

    for pdb_id in PDB_IDS:
        dry_rec = extract_receptor(structures[pdb_id], pdb_id)
        dry_receptors[pdb_id] = dry_rec

        lig_name = LIGAND_OF[pdb_id]
        lpdb = extract_aligned_ligand(structures[pdb_id], pdb_id, lig_name)
        lig_pdbs[pdb_id] = lpdb

        if lpdb:
            lsdf = extract_ligand_sdf(lpdb, pdb_id, lig_name)
            lig_sdfs[pdb_id] = lsdf

        waters = get_pharmacophore_waters(structures[pdb_id], pdb_id)
        wet_rec = create_wet_receptor(dry_rec, waters, pdb_id)
        wet_receptors[pdb_id] = wet_rec
        water_counts[pdb_id] = len(waters)

        sources = {}
        for w in waters:
            src = w["source"]
            sources[src] = sources.get(src, 0) + 1
        src_str = ", ".join(f"{k}={v}" for k, v in sorted(sources.items()))
        print(f"  {pdb_id}: {len(waters)} waters embedded ({src_str})")

    docking_pairs = []
    for rec_id in PDB_IDS:
        for lig_id in PDB_IDS:
            lig_name = LIGAND_OF[lig_id]
            is_self = (rec_id == lig_id)
            pair_type = "self" if is_self else "cross"
            docking_pairs.append((rec_id, lig_id, lig_name, pair_type))

    print(f"\n{'=' * 70}")
    print(f"Step 3: Running {len(docking_pairs)} docking pairs "
          f"(dry and wet)...")
    print(f"{'=' * 70}")

    all_results = []

    for rec_id, lig_id, lig_name, pair_type in docking_pairs:
        lig_sdf = lig_sdfs.get(lig_id)
        lig_pdb = lig_pdbs.get(lig_id)
        autobox_ref = lig_pdbs.get(rec_id)

        if not lig_sdf or not lig_pdb or not autobox_ref:
            print(f"\n  {lig_name} -> {rec_id}: SKIPPED (missing files)")
            continue

        print(f"\n  {lig_name} ({lig_id}) -> {rec_id} ({pair_type}):")

        result = {
            "receptor": rec_id,
            "ligand_source": lig_id,
            "ligand": lig_name,
            "pair_type": pair_type,
        }

        for mode in ["dry", "wet"]:
            rec_path = dry_receptors[rec_id] if mode == "dry" \
                else wet_receptors[rec_id]
            out_sdf = os.path.join(
                OUTPUT_DIR,
                f"{rec_id}_{lig_id}_{lig_name}_{mode}.sdf")

            success = run_gnina_dock(
                rec_path, lig_sdf, autobox_ref, out_sdf,
                label=f"{lig_name}->{rec_id} ({mode})")

            if success:
                rmsd = calc_rmsd(out_sdf, lig_pdb)
                cnn, vina = get_gnina_scores(out_sdf)

                include_w = (mode == "wet")
                ifp = compute_ifp(rec_path, out_sdf, rec_id,
                                  include_waters=include_w)

                crystal_sdf_path = os.path.join(
                    OUTPUT_DIR,
                    f"{rec_id}_{LIGAND_OF[rec_id]}_{LIGAND_OF[rec_id]}_"
                    f"{mode}.sdf")
                if rec_id == lig_id:
                    crystal_ifp = ifp.copy()
                    result[f"{mode}_crystal_ifp"] = ifp.tolist()

                rmsd_s = f"{rmsd:.3f}" if rmsd else "N/A"
                cnn_s = f"{cnn:.4f}" if cnn else "N/A"
                vina_s = f"{vina:.2f}" if vina else "N/A"
                bits_on = int(np.sum(ifp))

                print(f"    {mode.upper():>3s}: RMSD={rmsd_s}, CNN={cnn_s}, "
                      f"Vina={vina_s}, IFP bits={bits_on}")

                result[f"{mode}_rmsd"] = rmsd
                result[f"{mode}_cnn"] = cnn
                result[f"{mode}_vina"] = vina
                result[f"{mode}_ifp"] = ifp.tolist()
                result[f"{mode}_bits_on"] = bits_on
            else:
                print(f"    {mode.upper():>3s}: FAILED")
                result[f"{mode}_rmsd"] = None
                result[f"{mode}_cnn"] = None
                result[f"{mode}_vina"] = None

        all_results.append(result)

    print(f"\n{'=' * 70}")
    print("Step 4: Computing Tanimoto scores...")
    print(f"{'=' * 70}")

    crystal_ifps = {"dry": {}, "wet": {}}
    for r in all_results:
        if r["pair_type"] == "self":
            rec = r["receptor"]
            for mode in ["dry", "wet"]:
                ifp = r.get(f"{mode}_ifp")
                if ifp:
                    crystal_ifps[mode][rec] = np.array(ifp)

    for r in all_results:
        rec = r["receptor"]
        for mode in ["dry", "wet"]:
            ifp = r.get(f"{mode}_ifp")
            ref_ifp = crystal_ifps[mode].get(rec)
            if ifp is not None and ref_ifp is not None:
                tc = tanimoto(np.array(ifp), ref_ifp)
                r[f"{mode}_tc"] = tc
            else:
                r[f"{mode}_tc"] = None

    print(f"\n  {'Ligand':>8s} -> {'Receptor':>8s} {'Type':>6s}  "
          f"{'Dry Tc':>7s} {'Wet Tc':>7s} {'Delta':>7s}  "
          f"{'Dry RMSD':>9s} {'Wet RMSD':>9s}")
    print(f"  {'-' * 75}")

    for r in all_results:
        lig_label = f"{r['ligand']}({r['ligand_source'][:4]})"
        rec_label = r["receptor"]
        ptype = r["pair_type"]

        dry_tc = r.get("dry_tc")
        wet_tc = r.get("wet_tc")
        dry_rmsd = r.get("dry_rmsd")
        wet_rmsd = r.get("wet_rmsd")

        dry_tc_s = f"{dry_tc:.4f}" if dry_tc is not None else "  N/A"
        wet_tc_s = f"{wet_tc:.4f}" if wet_tc is not None else "  N/A"
        delta_s = f"{wet_tc - dry_tc:+.4f}" if (dry_tc is not None
                  and wet_tc is not None) else "  N/A"
        dry_r_s = f"{dry_rmsd:.3f}" if dry_rmsd else "  N/A"
        wet_r_s = f"{wet_rmsd:.3f}" if wet_rmsd else "  N/A"

        marker = ""
        if dry_tc is not None and wet_tc is not None:
            if wet_tc > dry_tc:
                marker = " (+)"
            elif wet_tc < dry_tc:
                marker = " (-)"

        print(f"  {lig_label:>8s} -> {rec_label:>8s} {ptype:>6s}  "
              f"{dry_tc_s:>7s} {wet_tc_s:>7s} {delta_s:>7s}{marker}  "
              f"{dry_r_s:>9s} {wet_r_s:>9s}")

    print(f"\n{'=' * 70}")
    print("Step 5: Summary")
    print(f"{'=' * 70}\n")

    self_results = [r for r in all_results if r["pair_type"] == "self"]
    cross_results = [r for r in all_results if r["pair_type"] == "cross"]

    for label, subset in [("Self-docking", self_results),
                          ("Cross-docking", cross_results)]:
        dry_tcs = [r["dry_tc"] for r in subset if r.get("dry_tc") is not None]
        wet_tcs = [r["wet_tc"] for r in subset if r.get("wet_tc") is not None]
        dry_rmsds = [r["dry_rmsd"] for r in subset
                     if r.get("dry_rmsd") is not None]
        wet_rmsds = [r["wet_rmsd"] for r in subset
                     if r.get("wet_rmsd") is not None]

        if dry_tcs and wet_tcs:
            print(f"  {label}:")
            print(f"    Dry  mean Tc: {np.mean(dry_tcs):.4f}, "
                  f"mean RMSD: {np.mean(dry_rmsds):.3f} A")
            print(f"    Wet  mean Tc: {np.mean(wet_tcs):.4f}, "
                  f"mean RMSD: {np.mean(wet_rmsds):.3f} A")
            delta_tc = np.mean(wet_tcs) - np.mean(dry_tcs)
            delta_rmsd = np.mean(wet_rmsds) - np.mean(dry_rmsds)
            print(f"    Delta Tc: {delta_tc:+.4f}, "
                  f"Delta RMSD: {delta_rmsd:+.3f} A")

            improved = sum(1 for d, w in zip(dry_tcs, wet_tcs) if w > d)
            print(f"    Improved: {improved}/{len(dry_tcs)}")
            print()

    all_dry_tc = [r["dry_tc"] for r in all_results
                  if r.get("dry_tc") is not None]
    all_wet_tc = [r["wet_tc"] for r in all_results
                  if r.get("wet_tc") is not None]

    if all_dry_tc and all_wet_tc:
        dry_pass = sum(1 for tc in all_dry_tc if tc >= TC_THRESHOLD)
        wet_pass = sum(1 for tc in all_wet_tc if tc >= TC_THRESHOLD)
        print(f"  Tanimoto >= {TC_THRESHOLD} threshold:")
        print(f"    Dry: {dry_pass}/{len(all_dry_tc)} pass")
        print(f"    Wet: {wet_pass}/{len(all_wet_tc)} pass")

        overall_delta = np.mean(all_wet_tc) - np.mean(all_dry_tc)
        verdict = "VALIDATED" if overall_delta > 0 else "NOT VALIDATED"
        print(f"\n  Overall delta Tc: {overall_delta:+.4f}")
        print(f"  Water placement: {verdict}")
        if overall_delta > 0:
            print(f"  Conserved waters IMPROVE interaction fidelity.")
            print(f"  Include them as explicit receptor waters in AutoGrow4.")
        else:
            print(f"  Waters did not improve scores. Consider flexible")
            print(f"  water treatment instead of fixed positions.")

    results_json = {
        "parameters": {
            "exhaustiveness": EXHAUSTIVENESS,
            "tc_threshold": TC_THRESHOLD,
            "ifp_types": IFP_TYPES,
        },
        "water_counts": water_counts,
        "docking_results": [],
    }
    for r in all_results:
        entry = {k: v for k, v in r.items()
                 if not k.endswith("_ifp") and not k.endswith("crystal_ifp")}
        results_json["docking_results"].append(entry)

    json_path = os.path.join(OUTPUT_DIR,
                             "water_validated_results.json")
    with open(json_path, "w") as f:
        json.dump(results_json, f, indent=2, default=str)
    print(f"\nJSON: {json_path}")

    report_path = os.path.join(OUTPUT_DIR,
                               "water_validated_report.txt")
    _write_report(report_path, all_results, water_counts)
    print(f"Report: {report_path}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        _generate_plots(all_results, water_counts)
    except ImportError:
        print("matplotlib not available, skipping plots")


def _write_report(path, all_results, water_counts):
    with open(path, "w") as f:
        f.write("WATER-VALIDATED CROSS-DOCKING - Task 6\n")
        f.write(f"{'=' * 60}\n\n")

        f.write("Embedded waters per receptor:\n")
        for pdb_id in PDB_IDS:
            f.write(f"  {pdb_id}: {water_counts.get(pdb_id, 0)} waters\n")

        f.write(f"\n{'=' * 60}\n")
        f.write(f"{'Ligand':>12s} -> {'Receptor':>8s} {'Type':>6s}  "
                f"{'Dry Tc':>7s} {'Wet Tc':>7s} {'Delta':>7s}  "
                f"{'Dry RMSD':>9s} {'Wet RMSD':>9s}\n")
        f.write(f"{'-' * 75}\n")

        for r in all_results:
            lig_label = f"{r['ligand']}({r['ligand_source'][:4]})"
            dry_tc = r.get("dry_tc")
            wet_tc = r.get("wet_tc")
            dry_rmsd = r.get("dry_rmsd")
            wet_rmsd = r.get("wet_rmsd")

            f.write(f"{lig_label:>12s} -> {r['receptor']:>8s} "
                    f"{r['pair_type']:>6s}  ")
            f.write(f"{dry_tc:.4f} " if dry_tc is not None else "  N/A  ")
            f.write(f"{wet_tc:.4f} " if wet_tc is not None else "  N/A  ")
            if dry_tc is not None and wet_tc is not None:
                f.write(f"{wet_tc - dry_tc:+.4f}  ")
            else:
                f.write("  N/A    ")
            f.write(f"{dry_rmsd:.3f} " if dry_rmsd else "  N/A  ")
            f.write(f"{wet_rmsd:.3f}\n" if wet_rmsd else "  N/A\n")

        self_r = [r for r in all_results if r["pair_type"] == "self"]
        cross_r = [r for r in all_results if r["pair_type"] == "cross"]

        f.write(f"\n{'=' * 60}\n")
        f.write("Summary:\n\n")

        for label, subset in [("Self-docking", self_r),
                              ("Cross-docking", cross_r)]:
            dry_tcs = [r["dry_tc"] for r in subset
                       if r.get("dry_tc") is not None]
            wet_tcs = [r["wet_tc"] for r in subset
                       if r.get("wet_tc") is not None]
            if dry_tcs and wet_tcs:
                f.write(f"  {label}:\n")
                f.write(f"    Dry mean Tc: {np.mean(dry_tcs):.4f}\n")
                f.write(f"    Wet mean Tc: {np.mean(wet_tcs):.4f}\n")
                f.write(f"    Delta: {np.mean(wet_tcs) - np.mean(dry_tcs):+.4f}\n\n")


def _generate_plots(all_results, water_counts):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    fig.suptitle("Water-Validated Cross-Docking - Task 6\n"
                 "Dry vs Wet Receptor Comparison",
                 fontsize=14, fontweight="bold")

    labels = []
    dry_tcs = []
    wet_tcs = []
    pair_types = []

    for r in all_results:
        dry_tc = r.get("dry_tc")
        wet_tc = r.get("wet_tc")
        if dry_tc is not None and wet_tc is not None:
            lig_label = f"{r['ligand']}\n({r['ligand_source'][:4]})\n->{r['receptor']}"
            labels.append(lig_label)
            dry_tcs.append(dry_tc)
            wet_tcs.append(wet_tc)
            pair_types.append(r["pair_type"])

    if labels:
        x = np.arange(len(labels))
        width = 0.35

        ax = axes[0]
        colors_dry = ["#90CAF9" if t == "self" else "#CE93D8"
                      for t in pair_types]
        colors_wet = ["#2196F3" if t == "self" else "#9C27B0"
                      for t in pair_types]
        ax.bar(x - width / 2, dry_tcs, width, label="Dry",
               color=colors_dry, edgecolor="black", linewidth=0.5)
        ax.bar(x + width / 2, wet_tcs, width, label="Wet",
               color=colors_wet, edgecolor="black", linewidth=0.5)
        ax.axhline(y=TC_THRESHOLD, color="red", linestyle="--",
                   linewidth=1.5, label=f"Tc={TC_THRESHOLD}")
        ax.set_ylabel("Tanimoto Coefficient", fontsize=11)
        ax.set_title("IFP Tanimoto: Dry vs Wet", fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7)
        ax.legend(fontsize=9)
        ax.grid(axis="y", alpha=0.3)
        ax.set_ylim(0, 1.1)

    deltas = [w - d for d, w in zip(dry_tcs, wet_tcs)]
    if deltas:
        ax2 = axes[1]
        colors = ["#4CAF50" if d > 0 else "#F44336" for d in deltas]
        ax2.bar(range(len(deltas)), deltas, color=colors,
                edgecolor="black", linewidth=0.5)
        ax2.axhline(y=0, color="black", linewidth=1)
        ax2.set_ylabel("Delta Tanimoto (Wet - Dry)", fontsize=11)
        ax2.set_title("Water Impact on IFP", fontweight="bold")
        ax2.set_xticks(range(len(labels)))
        ax2.set_xticklabels(labels, fontsize=7)
        ax2.grid(axis="y", alpha=0.3)

        for i, d in enumerate(deltas):
            ax2.text(i, d + 0.002 * (1 if d >= 0 else -1),
                     f"{d:+.3f}", ha="center", fontsize=7,
                     fontweight="bold")

    dry_rmsds = []
    wet_rmsds = []
    rmsd_labels = []
    for r in all_results:
        dr = r.get("dry_rmsd")
        wr = r.get("wet_rmsd")
        if dr is not None and wr is not None:
            dry_rmsds.append(dr)
            wet_rmsds.append(wr)
            rmsd_labels.append(
                f"{r['ligand']}\n({r['ligand_source'][:4]})\n->{r['receptor']}")

    if rmsd_labels:
        ax3 = axes[2]
        x = np.arange(len(rmsd_labels))
        ax3.bar(x - width / 2, dry_rmsds, width, label="Dry",
                color="#FFCC80", edgecolor="black", linewidth=0.5)
        ax3.bar(x + width / 2, wet_rmsds, width, label="Wet",
                color="#FF9800", edgecolor="black", linewidth=0.5)
        ax3.axhline(y=0.5, color="red", linestyle="--",
                    linewidth=1.5, label="0.5 A")
        ax3.set_ylabel("RMSD (A)", fontsize=11)
        ax3.set_title("RMSD: Dry vs Wet", fontweight="bold")
        ax3.set_xticks(x)
        ax3.set_xticklabels(rmsd_labels, fontsize=7)
        ax3.legend(fontsize=9)
        ax3.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    save_path = os.path.join(OUTPUT_DIR, "water_validated_comparison.png")
    plt.savefig(save_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Comparison plot: {save_path}")


if __name__ == "__main__":
    main()
