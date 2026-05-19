"""
GNINA Redocking Pose Visualization
Paste this into a Jupyter notebook cell on Expanse.

Shows top GNINA docked pose vs crystal ligand for all 4 conditions:
  - 2ZSH GA3 pH 7.2/7.4
  - 2ZSI GA4 pH 7.2/7.4

Active site residues displayed as functional group sticks:
  119, 24, 27, 31, 126, 323, 239, 319, 116, 191, 244, 127, 238, 218

Color key:
  Green  = Crystal ligand (reference position from PDB)
  Magenta = GNINA top docked pose
  Sticks  = Active site residue side chains (element coloring)
  White   = Protein backbone (transparent cartoon)
"""

import json, os, time
import nglview as nv
from IPython.display import display, HTML

ACTIVE_SITE = [24, 27, 31, 116, 119, 126, 127, 191, 218, 238, 239, 244, 319, 323]
BASE = os.path.expanduser("~/final_project/redocking_gnina")
RESULTS_JSON = os.path.join(BASE, "redocking_results.json")

with open(RESULTS_JSON) as f:
    all_results = json.load(f)

crystal_conditions = [c for c in all_results if c.get("ligand_source") == "crystal"]

display(HTML("""
<div style='background:#1a1a2e; color:#e0e0e0; padding:15px; border-radius:8px; margin-bottom:15px;'>
    <h2 style='margin:0 0 10px 0; color:#00ff88;'>GNINA Redocking Pose Validation</h2>
    <p style='margin:5px 0;'>Comparing top GNINA docked pose against crystallographic ligand position.</p>
    <p style='margin:5px 0;'>
        <span style='color:#00ff88; font-size:1.3em;'>&#9632;</span> Crystal ligand (reference) &nbsp;&nbsp;
        <span style='color:#ff44aa; font-size:1.3em;'>&#9632;</span> GNINA docked pose &nbsp;&nbsp;
        <span style='color:#88ccff; font-size:1.3em;'>&#9632;</span> Active site residues (element coloring)
    </p>
    <p style='margin:5px 0; font-size:0.9em; color:#aaa;'>
        Active site: """ + ", ".join(str(r) for r in ACTIVE_SITE) + """
    </p>
</div>
"""))

print(f"Loaded {len(crystal_conditions)} conditions from {RESULTS_JSON}\n")

for cond in crystal_conditions:
    pdb_id = cond["pdb_id"]
    lig_name = cond["ligand_name"]
    ph = cond["ph"]

    sweep = cond.get("sweep_results", [])
    if not sweep:
        print(f"  No sweep results for {pdb_id} pH {ph}, skipping\n")
        continue

    best = min(sweep, key=lambda r: r.get("best_rmsd", 999))
    docked_file = best.get("docked_file", "")
    rmsd = best.get("best_rmsd", 999)
    vina = best.get("best_vina_score")
    cnn = best.get("best_cnn_score")

    pdb_file = os.path.join(BASE, pdb_id, f"{pdb_id}.pdb")
    crystal_lig = os.path.join(BASE, pdb_id,
                               f"{pdb_id}_ligand_{lig_name}.pdb")

    missing = False
    for fp, label in [(pdb_file, "Original PDB"),
                      (crystal_lig, "Crystal ligand"),
                      (docked_file, "Docked SDF")]:
        if not fp or not os.path.exists(fp):
            print(f"  MISSING {label}: {fp}")
            missing = True
    if missing:
        print()
        continue

    with open(pdb_file) as f:
        pdb_text = f.read()
    with open(crystal_lig) as f:
        crystal_text = f.read()
    with open(docked_file) as f:
        sdf_all = f.read()

    parts = sdf_all.split("$$$$")
    first_pose = parts[0] + "\n$$$$\n" if parts[0].strip() else ""
    if not first_pose:
        print(f"  Empty docked SDF for {pdb_id} pH {ph}\n")
        continue

    rmsd_s = f"{rmsd:.3f}" if isinstance(rmsd, float) else str(rmsd)
    vina_s = f"{vina:.2f} kcal/mol" if vina is not None else "N/A"
    cnn_s = f"{cnn:.4f}" if cnn is not None else "N/A"

    if best.get("method") == "autobox":
        method_s = f"autobox padding = {best['padding']:.0f} A"
    else:
        method_s = f"manual box = {best.get('size_x', '?')} A"

    display(HTML(f"""
    <div style='background:#0d1117; color:#e0e0e0; padding:12px; margin:15px 0 5px 0;
                border-radius:8px; border-left:4px solid #00ff88;'>
        <h3 style='margin:0; color:#00ff88;'>{pdb_id} | {lig_name} | pH {ph}</h3>
        <table style='margin:8px 0; color:#e0e0e0;'>
            <tr><td style='padding:2px 15px 2px 0;'>RMSD vs Crystal:</td>
                <td style='color:#ffcc00; font-weight:bold;'>{rmsd_s} A</td></tr>
            <tr><td style='padding:2px 15px 2px 0;'>Vina Score:</td>
                <td>{vina_s}</td></tr>
            <tr><td style='padding:2px 15px 2px 0;'>CNN Score:</td>
                <td>{cnn_s}</td></tr>
            <tr><td style='padding:2px 15px 2px 0;'>Best Config:</td>
                <td>{method_s}</td></tr>
        </table>
    </div>
    """))

    view = nv.NGLWidget()
    view._remote_call("setSize", target="Widget", args=["100%", "550px"])

    # Component 0: Full PDB structure (receptor)
    view.add_component(nv.TextStructure(pdb_text, ext="pdb"),
                       default_representation=False, name="receptor")
    time.sleep(0.3)

    view.add_cartoon(selection="protein", color="white", opacity=0.12,
                     component=0)

    res_sel = ("(" + " or ".join(str(r) for r in ACTIVE_SITE)
               + ") and sidechainAttached")
    view.add_licorice(selection=res_sel, color="element", radius=0.18,
                      component=0)

    label_sel = ("(" + " or ".join(str(r) for r in ACTIVE_SITE)
                 + ") and .CA")
    view.add_label(selection=label_sel, color="white",
                   labelType="format",
                   labelFormat="%(resname)s%(resno)s",
                   labelGrouping="residue", component=0,
                   zOffset=2.0, fixedSize=True,
                   attachment="middle_center",
                   showBackground=True,
                   backgroundColor="black",
                   backgroundOpacity=0.6)

    # Component 1: Crystal ligand (green)
    view.add_component(nv.TextStructure(crystal_text, ext="pdb"),
                       default_representation=False, name="crystal_ligand")
    time.sleep(0.3)
    view.add_ball_and_stick(selection="all", color="green",
                            aspectRatio=2.5, component=1)

    # Component 2: GNINA docked pose (magenta)
    view.add_component(nv.TextStructure(first_pose, ext="sdf"),
                       default_representation=False, name="gnina_docked")
    time.sleep(0.3)
    view.add_ball_and_stick(selection="all", color="#ff44aa",
                            aspectRatio=2.5, component=2)

    view.center(selection="all", component=1)

    display(view)
    print(f"  Rotate to inspect overlap. Green = crystal, Magenta = docked.\n")

display(HTML("""
<div style='background:#1a1a2e; color:#aaa; padding:10px; border-radius:8px; margin-top:15px;
            font-size:0.9em;'>
    <b>Interpretation:</b> If green and magenta overlap closely, GNINA successfully
    reproduced the crystallographic binding mode. Sub-0.5 A RMSD means the poses
    are virtually identical - the CNN scoring correctly penalizes flipped orientations
    that Vina scored identically.
</div>
"""))
