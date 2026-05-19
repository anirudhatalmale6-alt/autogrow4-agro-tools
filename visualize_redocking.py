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
<pre style='background:#ffffff !important; color:#000000 !important; padding:15px; border-radius:8px;
            margin-bottom:15px; border:2px solid #228B22; font-family:monospace; font-size:14px;
            line-height:1.6; white-space:pre;'><span style='color:#006400 !important; font-size:22px; font-weight:bold;'>GNINA Redocking Pose Validation</span>

Comparing top GNINA docked pose against crystallographic ligand position.

  <span style='color:#228B22 !important; font-size:16px;'>&#9632;</span> Crystal ligand (reference)
  <span style='color:#ff44aa !important; font-size:16px;'>&#9632;</span> GNINA docked pose
  <span style='color:#888 !important; font-size:16px;'>&#9632;</span> Active site residues (element coloring)

Active site: """ + ", ".join(str(r) for r in ACTIVE_SITE) + """</pre>
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
    vina_s = f"{vina:.2f}" if vina is not None else "N/A"
    cnn_s = f"{cnn:.4f}" if cnn is not None else "N/A"
    center = cond.get("center", [0, 0, 0])
    cx, cy, cz = center[0], center[1], center[2]

    if best.get("method") == "autobox":
        sx = best.get("size_x")
        if sx is not None:
            dim_s = f"{sx:.1f} x {best['size_y']:.1f} x {best['size_z']:.1f}"
        else:
            p = best.get("padding", 4)
            dim_s = f"autobox (ligand + {p:.0f} A padding each side)"
        method_s = f"Autobox, padding = {best['padding']:.0f} A"
    else:
        sx = best.get("size_x", 0)
        sy = best.get("size_y", 0)
        sz = best.get("size_z", 0)
        dim_s = f"{sx:.1f} x {sy:.1f} x {sz:.1f} A"
        method_s = f"Manual box"

    sweep_rows = ""
    for r in sorted(sweep, key=lambda x: x.get("best_rmsd", 999)):
        r_rmsd = r.get("best_rmsd", 999)
        if r_rmsd >= 999:
            continue
        r_rmsd_s = f"{r_rmsd:.3f}"
        r_vina = r.get("best_vina_score")
        r_vina_s = f"{r_vina:.2f}" if r_vina is not None else "N/A"
        is_best = " BEST" if r is best else ""
        if r.get("method") == "autobox":
            r_dim = f"Autobox pad {r['padding']:.0f}A"
        else:
            r_dim = f"{r.get('size_x', 0):.0f}x{r.get('size_y', 0):.0f}x{r.get('size_z', 0):.0f} A"
        bg = "background:#d4edda !important;" if r is best else "background:#ffffff !important;"
        weight = "font-weight:bold;" if r is best else ""
        sweep_rows += (f"<tr style='{bg} {weight}'>"
                       f"<td style='padding:3px 12px; color:#000 !important; background:inherit;'>{r_dim}</td>"
                       f"<td style='padding:3px 12px; color:#000 !important; background:inherit;'>{r_rmsd_s}</td>"
                       f"<td style='padding:3px 12px; color:#000 !important; background:inherit;'>{r_vina_s}</td>"
                       f"<td style='padding:3px 12px; color:#006400 !important; background:inherit; font-weight:bold;'>{is_best}</td></tr>")

    info_text = f"""
    <pre style='background:#ffffff !important; color:#000000 !important; padding:16px; margin:15px 0 5px 0;
                border-radius:8px; border:2px solid #228B22; border-left:5px solid #228B22;
                font-family:monospace; font-size:14px; line-height:1.8; white-space:pre;'><span style='color:#006400 !important; font-size:20px; font-weight:bold;'>{pdb_id} | {lig_name} | pH {ph}</span>

  RMSD vs Crystal:   <span style='color:#cc0000 !important; font-weight:bold; font-size:16px;'>{rmsd_s} A</span>
  Vina Score:        {vina_s} kcal/mol
  CNN Score:         {cnn_s}
  Best Method:       {method_s}
  Box Dimensions:    {dim_s}
  Box Center:        ({cx:.2f}, {cy:.2f}, {cz:.2f})</pre>"""

    display(HTML(info_text))

    display(HTML(f"""
    <details style='margin:0 0 5px 0;'>
        <summary style='cursor:pointer; color:#0066cc !important; font-size:14px; font-weight:bold;
                        padding:6px; background:#f8f8f8 !important; border:1px solid #ccc; border-radius:4px;'>
            All configurations tested (click to expand)</summary>
        <table style='margin:8px 0; font-size:13px; border-collapse:collapse;
                      background:#ffffff !important; border:1px solid #ccc;'>
            <tr style='background:#e8e8e8 !important; border-bottom:2px solid #999;'>
                <th style='padding:5px 12px; text-align:left; color:#000 !important; background:#e8e8e8 !important;'>Configuration</th>
                <th style='padding:5px 12px; text-align:left; color:#000 !important; background:#e8e8e8 !important;'>RMSD (A)</th>
                <th style='padding:5px 12px; text-align:left; color:#000 !important; background:#e8e8e8 !important;'>Vina (kcal/mol)</th>
                <th style='padding:5px 12px; color:#000 !important; background:#e8e8e8 !important;'></th></tr>
            {sweep_rows}
        </table>
    </details>
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
<pre style='background:#ffffff !important; color:#000000 !important; padding:12px; border-radius:8px;
            margin-top:15px; font-size:13px; border:1px solid #228B22; font-family:monospace;
            line-height:1.5; white-space:pre-wrap;'>Interpretation: If green and magenta overlap closely, GNINA successfully reproduced the crystallographic binding mode. Sub-0.5 A RMSD means the poses are virtually identical - the CNN scoring correctly penalizes flipped orientations that Vina scored identically.</pre>
"""))
