"""
Analysis 06: Ligand H-Bond Frequency
Mirrors Figure 10 of Gődény et al. 2025.

H-bond criteria (same as paper):
  donor-acceptor distance <= 4.0 Å
  donor-H-acceptor angle  >= 120°

Outputs:
  results/ligand_hbonds_frequency.png   → Figure 10
  results/ligand_hbonds_timeseries.png
  results/ligand_hbonds.csv
"""

import os, glob
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TOP      = os.path.join(SCRIPT_DIR, "../../amoeba_ready_pdbs/H0_amoeba_ready.pdb")
TRAJ_DIR = os.path.join(SCRIPT_DIR, "../traj")
OUT      = os.path.join(SCRIPT_DIR, "results/06_ligand_hbonds")
os.makedirs(OUT, exist_ok=True)

TRAJS = sorted(glob.glob(os.path.join(TRAJ_DIR, "npt_*.dcd")), key=lambda p: int(os.path.basename(p).split("_")[-1].split(".")[0]))
if not TRAJS:
    raise FileNotFoundError(f"No npt_*.dcd found in {TRAJ_DIR}")

u = mda.Universe(TOP, TRAJS)
print(f"Ligand atoms (resname 308): {len(u.select_atoms('resname 308'))}")

# Ligand as donor → protein as acceptor
hba = HydrogenBondAnalysis(
    universe=u,
    donors_sel="resname 308",
    hydrogens_sel="resname 308 and name H*",
    acceptors_sel="protein and name N* O* S*",
    d_a_cutoff=4.0,
    d_h_a_angle_cutoff=120.0,
    update_selections=False,
)
hba.run(verbose=True)

# Protein as donor → ligand as acceptor
hba2 = HydrogenBondAnalysis(
    universe=u,
    donors_sel="protein",
    hydrogens_sel="protein and name H*",
    acceptors_sel="resname 308 and name N* O* S*",
    d_a_cutoff=4.0,
    d_h_a_angle_cutoff=120.0,
    update_selections=False,
)
hba2.run(verbose=True)

total_frames = u.trajectory.n_frames
resid_counts = defaultdict(int)

def parse_hbonds(hba_result, protein_role="acceptor"):
    if hba_result.results.hbonds is None or len(hba_result.results.hbonds) == 0:
        return
    for hb in hba_result.results.hbonds:
        frame, donor_idx, h_idx, acc_idx, dist, angle = hb
        atom = u.atoms[int(acc_idx if protein_role == "acceptor" else donor_idx)]
        resid_counts[f"{atom.resname}{atom.resid}"] += 1

parse_hbonds(hba, "acceptor")
parse_hbonds(hba2, "donor")

resid_freq  = {k: v / total_frames for k, v in resid_counts.items()}
sorted_res  = sorted(resid_freq.items(), key=lambda x: -x[1])

# H-bond count per frame
hb_per_frame = np.zeros(total_frames)
for arr in [hba.results.hbonds, hba2.results.hbonds]:
    if arr is not None and len(arr) > 0:
        hb_per_frame += np.bincount(arr[:, 0].astype(int), minlength=total_frames)

with open(os.path.join(OUT, "ligand_hbonds.csv"), "w") as f:
    f.write("residue,hbond_frequency\n")
    for res, freq in sorted_res:
        f.write(f"{res},{freq:.4f}\n")

times_ns = np.arange(total_frames) * (u.trajectory.dt / 1000)

# ── Plot 1: frequency bar (Figure 10) ────────────────────────────────────────
if sorted_res:
    labels = [r[0] for r in sorted_res[:20]]
    freqs  = [r[1] for r in sorted_res[:20]]
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(labels, freqs, color="steelblue", edgecolor="k", lw=0.5)
    ax.set_xlabel("Protein Residue")
    ax.set_ylabel("H-bond Frequency (fraction of frames)")
    ax.set_title("Ligand (308/amantadine) — Protein H-bond Frequency\n"
                 "AMOEBA Dipole Production  |  Figure 10 equivalent")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "ligand_hbonds_frequency.png"), dpi=150)
    plt.close()
    print("Saved: ligand_hbonds_frequency.png")
    print("Top 5 contacts:")
    for r, f in sorted_res[:5]:
        print(f"  {r}: {f*100:.1f}% of frames")
else:
    print("No H-bonds detected between ligand and protein.")

# ── Plot 2: timeseries ────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(times_ns, hb_per_frame, lw=0.7, color="darkorange")
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Number of H-bonds")
ax.set_title("Ligand–Protein H-bonds per Frame — AMOEBA Dipole Production")
ax.set_xlim(0, times_ns[-1])
fig.tight_layout()
fig.savefig(os.path.join(OUT, "ligand_hbonds_timeseries.png"), dpi=150)
plt.close()
print("Saved: ligand_hbonds_timeseries.png")
print("Done.")
