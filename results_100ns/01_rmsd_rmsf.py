"""
Analysis 01: RMSD and RMSF
Mirrors Figure 5 of Gődény et al. 2025 (Frontiers in Pharmacology).

Outputs:
  results/01_rmsd_rmsf/rmsd_backbone.png
  results/01_rmsd_rmsf/rmsf_per_residue.png
  results/01_rmsd_rmsf/rmsd.csv
  results/01_rmsd_rmsf/rmsf.csv
"""

import os, glob
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD, RMSF

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TOP    = os.path.join(SCRIPT_DIR, "../../amoeba_ready_pdbs/H0_amoeba_ready.pdb")
TRAJ_DIR = os.path.join(SCRIPT_DIR, "../traj")
OUT    = os.path.join(SCRIPT_DIR, "results/01_rmsd_rmsf")
os.makedirs(OUT, exist_ok=True)

# Auto-detect production trajectories (npt_*.dcd, sorted)
TRAJS = sorted(glob.glob(os.path.join(TRAJ_DIR, "npt_*.dcd")), key=lambda p: int(os.path.basename(p).split("_")[-1].split(".")[0]))
if not TRAJS:
    raise FileNotFoundError(f"No npt_*.dcd files found in {TRAJ_DIR}")
print(f"Using {len(TRAJS)} trajectory file(s): {[os.path.basename(t) for t in TRAJS]}")

u = mda.Universe(TOP, TRAJS)

# ── Fix periodic-boundary imaging artifact ────────────────────────────────────
# Near the end of npt_10.dcd, chain PROA's Ca atoms intermittently wrap across
# the periodic y-boundary (observed alternating between ~50.5 A and ~-0.3 A,
# exactly one box length apart), producing spurious ~20 A spikes in RMSD. This
# topology carries no protein backbone bonds (CONECT records only cover POPC
# lipid and ligand), so molecule/fragment-based wrapping isn't available;
# instead we remove box-length jumps by frame-to-frame continuity for all
# protein atoms before any RMSD/RMSF is computed.
_protein_all = u.select_atoms("protein")
_unwrap_state = {"prev": None, "shift": None}

def _unwrap_protein_jumps(ts):
    box = ts.dimensions[:3]
    pos = _protein_all.positions.copy()
    if _unwrap_state["prev"] is None:
        _unwrap_state["shift"] = np.zeros_like(pos)
    else:
        delta = pos - _unwrap_state["prev"]
        shift = _unwrap_state["shift"]
        shift = shift + np.where(delta < -box / 2, box, 0.0) \
                       + np.where(delta >  box / 2, -box, 0.0)
        _unwrap_state["shift"] = shift
    _unwrap_state["prev"] = pos
    _protein_all.positions = pos + _unwrap_state["shift"]
    return ts

u.trajectory.add_transformations(_unwrap_protein_jumps)

backbone = u.select_atoms("backbone")
ca       = u.select_atoms("name CA")

print(f"Loaded: {u.trajectory.n_frames} frames, {len(ca)} Cα atoms")

# ── RMSD ─────────────────────────────────────────────────────────────────────
rmsd_analysis = RMSD(backbone, backbone, select="backbone", ref_frame=0)
rmsd_analysis.run(verbose=True)

times     = rmsd_analysis.results.rmsd[:, 1] / 1000   # ps → ns
rmsd_vals = rmsd_analysis.results.rmsd[:, 2]           # Å

np.savetxt(os.path.join(OUT, "rmsd.csv"),
           np.column_stack([times, rmsd_vals]),
           header="time_ns,rmsd_A", delimiter=",", comments="")

mean_last_half = np.mean(rmsd_vals[len(rmsd_vals)//2:])
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(times, rmsd_vals, lw=0.8, color="steelblue")
ax.axhline(mean_last_half, color="red", ls="--", lw=1.2,
           label=f"Mean (last half): {mean_last_half:.2f} Å")
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Backbone RMSD (Å)")
ax.set_title("Backbone RMSD — AMOEBA Dipole Production")
ax.legend()
ax.set_xlim(0, times[-1])
fig.tight_layout()
fig.savefig(os.path.join(OUT, "rmsd_backbone.png"), dpi=150)
plt.close()
print("Saved: rmsd_backbone.png")

# ── RMSF ─────────────────────────────────────────────────────────────────────
rmsf_analysis = RMSF(ca).run(verbose=True)
rmsf_vals = rmsf_analysis.results.rmsf
resids    = ca.resids
resnames  = ca.resnames

np.savetxt(os.path.join(OUT, "rmsf.csv"),
           np.column_stack([resids, rmsf_vals]),
           header="resid,rmsf_A", delimiter=",", comments="",
           fmt=["%d", "%.4f"])

fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(resids, rmsf_vals, lw=0.9, color="darkorange")
his_mask = resnames == "HSD"
ax.scatter(resids[his_mask], rmsf_vals[his_mask],
           color="red", zorder=5, s=40, label="His37 (HSD)")
ax.set_xlabel("Residue ID")
ax.set_ylabel("Cα RMSF (Å)")
ax.set_title("Per-Residue RMSF — AMOEBA Dipole Production")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(OUT, "rmsf_per_residue.png"), dpi=150)
plt.close()
print("Saved: rmsf_per_residue.png")
print("Done.")
