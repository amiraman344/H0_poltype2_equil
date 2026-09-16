"""
Analysis 04: RDF g000 — Water around His37 Nitrogens
Mirrors Figure 8 (left, g000 panels) of Gődény et al. 2025.

Computes g(r) between His37 Nδ1/Nε2 and water O and H atoms.
For H0 system: all His37 are neutral (HSD), so:
  - Nδ1 carries the H (protonated δ-nitrogen)
  - Nε2 is deprotonated (lone pair)

This maps to the 3rd and 4th rows of Figure 8 left (Nδ and Nε panels).

Outputs:
  results/rdf_water_his_delta.png    → Figure 8, Nδ panels
  results/rdf_water_his_epsilon.png  → Figure 8, Nε panels
  results/rdf.csv
"""

import os, glob
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TOP      = os.path.join(SCRIPT_DIR, "../../amoeba_ready_pdbs/H0_amoeba_ready.pdb")
TRAJ_DIR = os.path.join(SCRIPT_DIR, "../traj")
OUT      = os.path.join(SCRIPT_DIR, "results/04_rdf_water_his")
os.makedirs(OUT, exist_ok=True)

TRAJS = sorted(glob.glob(os.path.join(TRAJ_DIR, "npt_*.dcd")), key=lambda p: int(os.path.basename(p).split("_")[-1].split(".")[0]))
if not TRAJS:
    raise FileNotFoundError(f"No npt_*.dcd found in {TRAJ_DIR}")

u = mda.Universe(TOP, TRAJS)

# His37 nitrogens (HSD = neutral, δ-protonated His)
his_nd1 = u.select_atoms("resname HID and name ND1")   # Nδ1 — carries H in HSD
his_ne2 = u.select_atoms("resname HID and name NE2")   # Nε2 — lone pair in HSD
wat_o   = u.select_atoms("resname HOH and name OH2")
wat_h   = u.select_atoms("resname HOH and name H1 H2")

print(f"His ND1: {len(his_nd1)}  His NE2: {len(his_ne2)}")
print(f"Water O: {len(wat_o)}  Water H: {len(wat_h)}")

nbins = 150
r_max = 8.0

def run_rdf(g1, g2, label):
    rdf = InterRDF(g1, g2, nbins=nbins, range=(0.5, r_max), verbose=True)
    rdf.run()
    return rdf.results.bins, rdf.results.rdf

print("\nRDF: water O around His Nδ1...")
bins_d, rdf_od = run_rdf(his_nd1, wat_o, "Nd1-O")
print("RDF: water O around His Nε2...")
bins_e, rdf_oe = run_rdf(his_ne2, wat_o, "Ne2-O")
print("RDF: water H around His Nδ1...")
_, rdf_hd = run_rdf(his_nd1, wat_h, "Nd1-H")
print("RDF: water H around His Nε2...")
_, rdf_he = run_rdf(his_ne2, wat_h, "Ne2-H")

np.savetxt(os.path.join(OUT, "rdf.csv"),
           np.column_stack([bins_d, rdf_od, rdf_oe, rdf_hd, rdf_he]),
           header="r_A,rdf_Nd1_O,rdf_Ne2_O,rdf_Nd1_H,rdf_Ne2_H",
           delimiter=",", comments="")
print("Saved: rdf.csv")

# ── Plot: Nδ1 ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(bins_d, rdf_od, lw=1.5, color="steelblue", label="g(r) Nδ1–O")
ax.plot(bins_d, rdf_hd, lw=1.5, color="darkorange", ls="--", label="g(r) Nδ1–H")
if len(rdf_od) > 0:
    peak = bins_d[np.argmax(rdf_od[:int(len(rdf_od)*0.5)])]
    ax.axvline(peak, color="steelblue", ls=":", lw=1, label=f"1st peak: {peak:.2f} Å")
ax.set_xlabel("r (Å)")
ax.set_ylabel("g(r)")
ax.set_title("RDF: Water around His37 Nδ1\n(H0 neutral His, AMOEBA)")
ax.legend()
ax.set_xlim(0.5, r_max)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "rdf_water_his_delta.png"), dpi=150)
plt.close()
print("Saved: rdf_water_his_delta.png")

# ── Plot: Nε2 ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(bins_e, rdf_oe, lw=1.5, color="crimson", label="g(r) Nε2–O")
ax.plot(bins_e, rdf_he, lw=1.5, color="purple", ls="--", label="g(r) Nε2–H")
if len(rdf_oe) > 0:
    peak = bins_e[np.argmax(rdf_oe[:int(len(rdf_oe)*0.5)])]
    ax.axvline(peak, color="crimson", ls=":", lw=1, label=f"1st peak: {peak:.2f} Å")
ax.set_xlabel("r (Å)")
ax.set_ylabel("g(r)")
ax.set_title("RDF: Water around His37 Nε2\n(H0 neutral His, AMOEBA)")
ax.legend()
ax.set_xlim(0.5, r_max)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "rdf_water_his_epsilon.png"), dpi=150)
plt.close()
print("Saved: rdf_water_his_epsilon.png")
print("Done.")
