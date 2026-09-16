"""
Analysis 08: Ligand MSD and Diffusion Coefficient
Mirrors Supplementary Figure S11 of Gődény et al. 2025.

Expected result: near-zero D → amantadine stays locked in binding site.

Outputs:
  results/ligand_msd.png
  results/ligand_diffusion.csv
"""

import os, glob
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from scipy.stats import linregress

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TOP      = os.path.join(SCRIPT_DIR, "../../amoeba_ready_pdbs/H0_amoeba_ready.pdb")
TRAJ_DIR = os.path.join(SCRIPT_DIR, "../traj")
OUT      = os.path.join(SCRIPT_DIR, "results/08_ligand_diffusion")
os.makedirs(OUT, exist_ok=True)

TRAJS = sorted(glob.glob(os.path.join(TRAJ_DIR, "npt_*.dcd")), key=lambda p: int(os.path.basename(p).split("_")[-1].split(".")[0]))
if not TRAJS:
    raise FileNotFoundError(f"No npt_*.dcd found in {TRAJ_DIR}")

u = mda.Universe(TOP, TRAJS)
ligand = u.select_atoms("resname 308")
print(f"Ligand: {len(ligand)} atoms")

print("Extracting ligand COM trajectory...")
com_traj = []
times_ps = []
for ts in u.trajectory:
    com_traj.append(ligand.center_of_mass())
    times_ps.append(ts.time)

com_traj = np.array(com_traj)
times_ns = np.array(times_ps) / 1000

# MSD via sliding window
N       = len(com_traj)
max_lag = N // 4
msd     = np.zeros(max_lag)
for lag in range(1, max_lag):
    disp      = com_traj[lag:] - com_traj[:-lag]
    msd[lag]  = np.mean(np.sum(disp**2, axis=1))

lag_times_ns = times_ns[:max_lag] - times_ns[0]

# Linear fit → D = slope / 6
fit_start = max_lag // 4
slope, intercept, r, _, _ = linregress(
    lag_times_ns[fit_start:] * 1e-9,
    msd[fit_start:] * 1e-20
)
D      = slope / 6.0
D_cm2s = D * 1e4
print(f"D = {D:.3e} m²/s = {D_cm2s:.3e} cm²/s  (expected near zero)")

np.savetxt(os.path.join(OUT, "ligand_diffusion.csv"),
           np.column_stack([lag_times_ns, msd]),
           header="lag_time_ns,msd_A2", delimiter=",", comments="")
with open(os.path.join(OUT, "ligand_diffusion.csv"), "a") as f:
    f.write(f"\n# D = {D:.4e} m2/s = {D_cm2s:.4e} cm2/s\n")
    f.write(f"# R2 = {r**2:.4f}\n")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].plot(lag_times_ns, msd, color="steelblue", lw=1.2, label="MSD (COM)")
fit_y = (slope * lag_times_ns[fit_start:] * 1e-9 + intercept) / 1e-20
axes[0].plot(lag_times_ns[fit_start:], fit_y,
             "r--", lw=1.5, label=f"D = {D:.2e} m²/s")
axes[0].set_xlabel("Lag time (ns)")
axes[0].set_ylabel("MSD (Å²)")
axes[0].set_title("Ligand (308) MSD")
axes[0].legend(fontsize=8)

axes[1].plot(times_ns, com_traj[:, 0] - com_traj[0, 0], lw=0.7, label="x", color="red")
axes[1].plot(times_ns, com_traj[:, 1] - com_traj[0, 1], lw=0.7, label="y", color="green")
axes[1].plot(times_ns, com_traj[:, 2] - com_traj[0, 2], lw=0.7, label="z", color="blue")
axes[1].set_xlabel("Time (ns)")
axes[1].set_ylabel("COM displacement (Å)")
axes[1].set_title("Ligand COM position")
axes[1].legend(fontsize=8)
axes[1].set_xlim(0, times_ns[-1])

fig.suptitle("Amantadine Diffusion — AMOEBA Dipole Production  |  Supp. Fig S11")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "ligand_msd.png"), dpi=150)
plt.close()
print("Saved: ligand_msd.png")
print("Done.")
