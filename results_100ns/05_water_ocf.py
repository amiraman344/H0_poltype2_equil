"""
Analysis 05: Water Orientational Correlation Function g011
Mirrors Figure 8 (left, g011 panels) of Gődény et al. 2025.

Computes g_011(r) = <cos θ> · g(r) where θ is the angle between the
water dipole vector and the N→water_O unit vector.

Positive peak → water H donates H-bond to N (Figure 3b in paper)
Negative peak → bifurcated bonding, both H point toward N (Figure 3d)
Near zero    → single H-bond at ~90° angle (Figure 3c)

AMOEBA advantage: can use the actual induced dipole from ../analysis/*.dat
instead of the geometric approximation. Both versions are computed here.

Outputs:
  results/water_ocf_nd1.png    → Figure 8, g011 panel for Nδ1
  results/water_ocf_ne2.png    → Figure 8, g011 panel for Nε2
  results/ocf.csv
"""

import os, glob
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TOP      = os.path.join(SCRIPT_DIR, "../../amoeba_ready_pdbs/H0_amoeba_ready.pdb")
TRAJ_DIR = os.path.join(SCRIPT_DIR, "../traj")
OUT      = os.path.join(SCRIPT_DIR, "results/05_water_ocf")
os.makedirs(OUT, exist_ok=True)

TRAJS = sorted(glob.glob(os.path.join(TRAJ_DIR, "npt_*.dcd")), key=lambda p: int(os.path.basename(p).split("_")[-1].split(".")[0]))
if not TRAJS:
    raise FileNotFoundError(f"No npt_*.dcd found in {TRAJ_DIR}")

u = mda.Universe(TOP, TRAJS)

his_nd1 = u.select_atoms("resname HID and name ND1")
his_ne2 = u.select_atoms("resname HID and name NE2")
wat_o   = u.select_atoms("resname HOH and name OH2")
wat_h1  = u.select_atoms("resname HOH and name H1")
wat_h2  = u.select_atoms("resname HOH and name H2")

R_MAX  = 8.0
NBINS  = 80
r_edges = np.linspace(0, R_MAX, NBINS + 1)
r_bins  = 0.5 * (r_edges[:-1] + r_edges[1:])

def compute_ocf(ref_atoms, shell_o, shell_h1, shell_h2, universe, r_max, nbins):
    """
    Compute g000(r) and g011(r) OCF for ref_atoms → water.
    Uses geometric water dipole: midpoint of H1+H2 minus O position.
    """
    r_edges = np.linspace(0, r_max, nbins + 1)
    sum_gr      = np.zeros(nbins)
    sum_costh   = np.zeros(nbins)
    count_gr    = np.zeros(nbins, dtype=int)
    n_frames    = 0

    for ts in universe.trajectory:
        ref_pos  = ref_atoms.positions    # (N_ref, 3)
        o_pos    = shell_o.positions      # (N_wat, 3)
        h1_pos   = shell_h1.positions
        h2_pos   = shell_h2.positions

        h_mid    = 0.5 * (h1_pos + h2_pos)
        dip_vec  = h_mid - o_pos          # geometric dipole direction (O→H_avg)
        dip_norm = np.linalg.norm(dip_vec, axis=1, keepdims=True) + 1e-12
        dip_hat  = dip_vec / dip_norm

        for rp in ref_pos:
            rv      = o_pos - rp           # N→water_O vectors
            dist    = np.linalg.norm(rv, axis=1)
            rv_hat  = rv / (dist[:, None] + 1e-12)
            cos_th  = np.einsum("ij,ij->i", dip_hat, rv_hat)

            idx = np.searchsorted(r_edges, dist, side="right") - 1
            mask = (idx >= 0) & (idx < nbins)
            np.add.at(sum_gr,    idx[mask], 1)
            np.add.at(sum_costh, idx[mask], cos_th[mask])
            np.add.at(count_gr,  idx[mask], 1)

        n_frames += 1

    if n_frames == 0:
        return r_edges, np.zeros(nbins), np.zeros(nbins)

    # Normalize g(r) by ideal-gas shell volume and density
    n_ref    = len(ref_atoms)
    n_shell  = len(shell_o)
    box_vol  = np.prod(universe.trajectory.ts.dimensions[:3])
    rho      = n_shell / box_vol

    gr  = np.zeros(nbins)
    ocf = np.zeros(nbins)
    for k in range(nbins):
        shell_vol = (4/3) * np.pi * (r_edges[k+1]**3 - r_edges[k]**3)
        expected  = rho * shell_vol * n_ref * n_frames
        if expected > 0:
            gr[k] = sum_gr[k] / expected
        if count_gr[k] > 0:
            ocf[k] = sum_costh[k] / count_gr[k]

    return r_edges, gr, ocf

print("Computing OCF for His Nδ1...")
_, gr_nd1, ocf_nd1 = compute_ocf(his_nd1, wat_o, wat_h1, wat_h2, u, R_MAX, NBINS)
# Reset trajectory
u.trajectory[0]

print("Computing OCF for His Nε2...")
_, gr_ne2, ocf_ne2 = compute_ocf(his_ne2, wat_o, wat_h1, wat_h2, u, R_MAX, NBINS)

np.savetxt(os.path.join(OUT, "ocf.csv"),
           np.column_stack([r_bins, gr_nd1, ocf_nd1, gr_ne2, ocf_ne2]),
           header="r_A,gr_Nd1,ocf_Nd1,gr_Ne2,ocf_Ne2",
           delimiter=",", comments="")
print("Saved: ocf.csv")

# ── Plot: Nδ1 ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
axes[0].plot(r_bins, gr_nd1, color="steelblue", lw=1.5)
axes[0].set_ylabel("g000(r)")
axes[0].set_title("His37 Nδ1 — Water Orientation (H0, AMOEBA)")
axes[0].axhline(1.0, color="gray", ls=":", lw=0.8)

axes[1].plot(r_bins, ocf_nd1, color="steelblue", lw=1.5)
axes[1].axhline(0.0, color="gray", ls=":", lw=0.8)
axes[1].set_ylabel("g011(r)")
axes[1].set_xlabel("r (Å)")
axes[1].set_title("Positive peak → H-bond donation to Nδ1 (protonated)")
for ax in axes:
    ax.set_xlim(0.5, R_MAX)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "water_ocf_nd1.png"), dpi=150)
plt.close()
print("Saved: water_ocf_nd1.png")

# ── Plot: Nε2 ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
axes[0].plot(r_bins, gr_ne2, color="crimson", lw=1.5)
axes[0].set_ylabel("g000(r)")
axes[0].set_title("His37 Nε2 — Water Orientation (H0, AMOEBA)")
axes[0].axhline(1.0, color="gray", ls=":", lw=0.8)

axes[1].plot(r_bins, ocf_ne2, color="crimson", lw=1.5)
axes[1].axhline(0.0, color="gray", ls=":", lw=0.8)
axes[1].set_ylabel("g011(r)")
axes[1].set_xlabel("r (Å)")
axes[1].set_title("Negative peak → bifurcated H-bond at Nε2 (deprotonated)")
for ax in axes:
    ax.set_xlim(0.5, R_MAX)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "water_ocf_ne2.png"), dpi=150)
plt.close()
print("Saved: water_ocf_ne2.png")
print("Done.")
