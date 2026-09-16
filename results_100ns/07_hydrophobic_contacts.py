"""
Analysis 07: Hydrophobic Contact Analysis — Ligand vs Protein
Mirrors Supplementary Figures S12-S15 of Gődény et al. 2025.

Cutoff: 4.5 Å heavy-atom distance between ligand and hydrophobic residues.

Outputs:
  results/hydrophobic_contacts_heatmap.png
  results/hydrophobic_contacts_timeseries.png
  results/hydrophobic_contacts.csv
"""

import os, glob
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TOP      = os.path.join(SCRIPT_DIR, "../../amoeba_ready_pdbs/H0_amoeba_ready.pdb")
TRAJ_DIR = os.path.join(SCRIPT_DIR, "../traj")
OUT      = os.path.join(SCRIPT_DIR, "results/07_hydrophobic_contacts")
os.makedirs(OUT, exist_ok=True)

TRAJS = sorted(glob.glob(os.path.join(TRAJ_DIR, "npt_*.dcd")), key=lambda p: int(os.path.basename(p).split("_")[-1].split(".")[0]))
if not TRAJS:
    raise FileNotFoundError(f"No npt_*.dcd found in {TRAJ_DIR}")

CUTOFF      = 4.5
HYDROPHOBIC = ("ALA", "VAL", "ILE", "LEU", "PRO", "TRP", "PHE", "MET", "GLY")

u = mda.Universe(TOP, TRAJS)
ligand   = u.select_atoms("resname 308 and not name H*")
hydropho = u.select_atoms(
    f"protein and not name H* and resname {' '.join(HYDROPHOBIC)}"
)
print(f"Ligand heavy atoms: {len(ligand)}")
print(f"Hydrophobic protein heavy atoms: {len(hydropho)}")

res_keys   = [f"{a.resname}{a.resid}" for a in hydropho]
unique_res = list(dict.fromkeys(res_keys))

contact_count  = defaultdict(int)
contacts_per_frame = []
total_frames = u.trajectory.n_frames

for ts in u.trajectory:
    dists = distance_array(ligand.positions, hydropho.positions, box=ts.dimensions)
    in_contact = dists < CUTOFF
    contacted  = set()
    for j, atom in enumerate(hydropho):
        if in_contact[:, j].any():
            key = f"{atom.resname}{atom.resid}"
            contact_count[key] += 1
            contacted.add(key)
    contacts_per_frame.append(len(contacted))
    if ts.frame % 200 == 0:
        print(f"  Frame {ts.frame}/{total_frames}")

contact_freq = {k: contact_count[k] / total_frames for k in unique_res}
sorted_res   = sorted(contact_freq.items(), key=lambda x: -x[1])

with open(os.path.join(OUT, "hydrophobic_contacts.csv"), "w") as f:
    f.write("residue,contact_frequency\n")
    for res, freq in sorted_res:
        f.write(f"{res},{freq:.4f}\n")

times_ns = np.arange(total_frames) * (u.trajectory.dt / 1000)

top20  = sorted_res[:20]
labels = [r[0] for r in top20]
freqs  = [r[1] for r in top20]
colors = ["tomato" if any(x in l for x in ("TRP","ALA","SER")) else "steelblue"
          for l in labels]

fig, ax = plt.subplots(figsize=(12, 5))
ax.bar(labels, freqs, color=colors, edgecolor="k", lw=0.5)
ax.set_xlabel("Protein Residue")
ax.set_ylabel("Contact Frequency (fraction of frames)")
ax.set_title("Ligand (308) Hydrophobic Contacts — AMOEBA Dipole Production\n"
             "Supp. Figs S12-S15 equivalent")
ax.tick_params(axis="x", rotation=45)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "hydrophobic_contacts_heatmap.png"), dpi=150)
plt.close()
print("Saved: hydrophobic_contacts_heatmap.png")

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(times_ns, contacts_per_frame, lw=0.7, color="seagreen")
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Residues in contact")
ax.set_title("Ligand Hydrophobic Contacts per Frame — AMOEBA Dipole Production")
ax.set_xlim(0, times_ns[-1])
fig.tight_layout()
fig.savefig(os.path.join(OUT, "hydrophobic_contacts_timeseries.png"), dpi=150)
plt.close()
print("Saved: hydrophobic_contacts_timeseries.png")
print("Done.")
