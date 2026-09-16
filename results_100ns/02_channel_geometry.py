"""
Analysis 02: Channel Cross-Distances and His37 Tetrad Geometry
Mirrors Figures 6 & 7 / 9 of Gődény et al. 2025.

BUG FIX vs original script: cross-distances are now computed as the
average of the two DIAGONAL pairs between OPPOSITE chains:
  d_top = ( dist(PROA_CA, PROC_CA) + dist(PROB_CA, PROD_CA) ) / 2
This matches the paper's definition (Figure 6: opposite Cα at top/bottom).
The old script used resids[:2] which picked adjacent atoms within one chain.

Chain layout (CHARMM-GUI SEGID naming):
  PROA  PROB
    \  /
    /  \
  PROD  PROC
  Diagonals: PROA↔PROC and PROB↔PROD

Residues:
  Top of channel (N-terminus): resid 22
  Bottom of channel (C-terminus): resid 46
  His37 gate: resid 37 (HSD, Nδ1 = ND1)

Outputs:
  results/channel_cross_distances.png   → Figure 6
  results/his37_tetrad_area.png         → Figure 7 (no ligand) / Figure 9 (with ligand)
  results/geometry.csv
"""

import os, glob
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TOP      = os.path.join(SCRIPT_DIR, "../../amoeba_ready_pdbs/H0_amoeba_ready.pdb")
TRAJ_DIR = os.path.join(SCRIPT_DIR, "../traj")
OUT      = os.path.join(SCRIPT_DIR, "results/02_channel_geometry")
os.makedirs(OUT, exist_ok=True)

TRAJS = sorted(glob.glob(os.path.join(TRAJ_DIR, "npt_*.dcd")), key=lambda p: int(os.path.basename(p).split("_")[-1].split(".")[0]))
if not TRAJS:
    raise FileNotFoundError(f"No npt_*.dcd found in {TRAJ_DIR}")
print(f"Trajectories: {[os.path.basename(t) for t in TRAJS]}")

u = mda.Universe(TOP, TRAJS)

# ── Fix periodic-boundary imaging artifact ────────────────────────────────────
# Near the end of npt_10.dcd, chain PROA's Ca atoms intermittently wrap across
# the periodic y-boundary (observed alternating between ~50.5 A and ~-0.3 A,
# exactly one box length apart), producing spurious spikes in the cross-
# distance and His37 tetrad geometry. This topology carries no protein
# backbone bonds (CONECT records only cover POPC lipid and ligand), so
# molecule/fragment-based wrapping isn't available; instead we remove
# box-length jumps by frame-to-frame continuity for all protein atoms before
# computing any distance.
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

# ── Chain segments (CHARMM-GUI SEGID convention) ──────────────────────────────
CHAINS = ["PROA", "PROB", "PROC", "PROD"]
TOP_RES = 22    # N-terminal residue id
BOT_RES = 46    # C-terminal residue id
HIS_RES = 37    # His37 gate residue id

# Verify chain CA counts
for seg in CHAINS:
    ca = u.select_atoms(f"segid {seg} and name CA")
    print(f"  {seg}: {len(ca)} CA atoms, resids {ca.resids[0]}–{ca.resids[-1]}")

# Select one CA per chain at top and bottom
ca_top = {seg: u.select_atoms(f"segid {seg} and name CA and resid {TOP_RES}")
          for seg in CHAINS}
ca_bot = {seg: u.select_atoms(f"segid {seg} and name CA and resid {BOT_RES}")
          for seg in CHAINS}

# His37 Nδ1 atoms (one per chain)
his_nd1 = u.select_atoms(f"resname HID and name ND1 and resid {HIS_RES}")
print(f"\nHSD ND1 atoms found: {len(his_nd1)}")
for a in his_nd1:
    print(f"  segid={a.segid} resid={a.resid} name={a.name}")

# ── Quadrilateral area from 4 3D points ──────────────────────────────────────
def quad_area(pts):
    """Area of a (possibly non-planar) quadrilateral via cross-product triangulation."""
    if len(pts) < 4:
        return np.nan
    p = pts - pts.mean(axis=0)
    _, _, vh = np.linalg.svd(p)
    proj = p @ vh[:2].T           # project onto best-fit plane
    x, y = proj[:, 0], proj[:, 1]
    n = len(x)
    return 0.5 * abs(sum(x[i]*y[(i+1)%n] - x[(i+1)%n]*y[i] for i in range(n)))

# ── Trajectory loop ───────────────────────────────────────────────────────────
times     = []
top_dists = []   # average diagonal cross-distance at channel top (Figure 6)
bot_dists = []   # average diagonal cross-distance at channel bottom (Figure 6)
his_areas = []   # area spanned by 4 Nδ1 atoms (Figure 7/9 left)
his_diag1 = []   # PROA↔PROC Nδ1 distance (Figure 7/9 right)
his_diag2 = []   # PROB↔PROD Nδ1 distance (Figure 7/9 right)

for ts in u.trajectory:
    times.append(ts.time / 1000)   # ps → ns

    # ── Cross-distances: average of two opposite-chain diagonal pairs ──────────
    # Paper: PROA↔PROC and PROB↔PROD diagonals (Figure 6)
    ok_top = all(len(ca_top[s]) == 1 for s in CHAINS)
    ok_bot = all(len(ca_bot[s]) == 1 for s in CHAINS)

    if ok_top:
        d_AC = np.linalg.norm(ca_top["PROA"].positions[0] - ca_top["PROC"].positions[0])
        d_BD = np.linalg.norm(ca_top["PROB"].positions[0] - ca_top["PROD"].positions[0])
        top_dists.append((d_AC + d_BD) / 2)
    else:
        top_dists.append(np.nan)

    if ok_bot:
        d_AC = np.linalg.norm(ca_bot["PROA"].positions[0] - ca_bot["PROC"].positions[0])
        d_BD = np.linalg.norm(ca_bot["PROB"].positions[0] - ca_bot["PROD"].positions[0])
        bot_dists.append((d_AC + d_BD) / 2)
    else:
        bot_dists.append(np.nan)

    # ── His37 tetrad geometry ─────────────────────────────────────────────────
    nd1_pos = his_nd1.positions
    if len(nd1_pos) == 4:
        # Order: PROA, PROB, PROC, PROD (sorted by segid → same as chain order)
        his_areas.append(quad_area(nd1_pos))
        # Diagonal 1: PROA(0)↔PROC(2), Diagonal 2: PROB(1)↔PROD(3)
        his_diag1.append(np.linalg.norm(nd1_pos[0] - nd1_pos[2]))
        his_diag2.append(np.linalg.norm(nd1_pos[1] - nd1_pos[3]))
    else:
        his_areas.append(np.nan)
        his_diag1.append(np.nan)
        his_diag2.append(np.nan)

times     = np.array(times)
top_dists = np.array(top_dists)
bot_dists = np.array(bot_dists)
his_areas = np.array(his_areas)
his_diag1 = np.array(his_diag1)
his_diag2 = np.array(his_diag2)

# ── Print averages for comparison with paper ─────────────────────────────────
print(f"\nMean cross-distance (top Cα): {np.nanmean(top_dists):.1f} Å  "
      f"(paper pol/pol H0 top: ~20.1 Å)")
print(f"Mean cross-distance (bot Cα): {np.nanmean(bot_dists):.1f} Å  "
      f"(paper pol/pol H0 bot: ~22.4 Å)")
print(f"Mean His37 Nδ area:           {np.nanmean(his_areas):.1f} Å²  "
      f"(paper pol/pol H0+lig: ~25.8 Å²)")
print(f"Mean His37 Nδ diag1 (A↔C):   {np.nanmean(his_diag1):.1f} Å")
print(f"Mean His37 Nδ diag2 (B↔D):   {np.nanmean(his_diag2):.1f} Å  "
      f"(paper pol/pol H0: ~7.1 Å)")

# ── Save CSV ──────────────────────────────────────────────────────────────────
np.savetxt(os.path.join(OUT, "geometry.csv"),
           np.column_stack([times, top_dists, bot_dists,
                            his_areas, his_diag1, his_diag2]),
           header="time_ns,top_dist_A,bot_dist_A,his_area_A2,his_diag1_A,his_diag2_A",
           delimiter=",", comments="")
print("\nSaved: geometry.csv")

# ── Plot 1: Cα cross-distances ────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(times, top_dists, label=f"Top Cα (res {TOP_RES}) avg diagonal",
        color="steelblue", lw=0.8)
ax.plot(times, bot_dists, label=f"Bottom Cα (res {BOT_RES}) avg diagonal",
        color="darkorange", lw=0.8)
ax.axhline(np.nanmean(top_dists), color="steelblue", ls="--", lw=1,
           alpha=0.7, label=f"Top mean: {np.nanmean(top_dists):.1f} Å")
ax.axhline(np.nanmean(bot_dists), color="darkorange", ls="--", lw=1,
           alpha=0.7, label=f"Bot mean: {np.nanmean(bot_dists):.1f} Å")
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Cross-distance (Å)")
ax.set_title("Channel Cα Diagonal Cross-Distances — AMOEBA Dipole Production\n"
             "(PROA↔PROC and PROB↔PROD averaged)")
ax.legend(fontsize=8)
ax.set_xlim(0, times[-1])
fig.tight_layout()
fig.savefig(os.path.join(OUT, "channel_cross_distances.png"), dpi=150)
plt.close()
print("Saved: channel_cross_distances.png")

# ── Plot 2: His37 tetrad ──────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].plot(times, his_areas, color="crimson", lw=0.8)
axes[0].axhline(np.nanmean(his_areas), color="k", ls="--", lw=1,
                label=f"Mean: {np.nanmean(his_areas):.1f} Å²")
axes[0].set_title(f"His37 Nδ1 Tetrad Area  (paper H0+lig pol/pol: ~25.8 Å²)")
axes[0].set_xlabel("Time (ns)")
axes[0].set_ylabel("Area (Å²)")
axes[0].legend()
axes[0].set_xlim(0, times[-1])

axes[1].plot(times, his_diag1, label="PROA↔PROC", color="purple", lw=0.8)
axes[1].plot(times, his_diag2, label="PROB↔PROD", color="teal", lw=0.8)
axes[1].axhline(np.nanmean(his_diag1), color="purple", ls="--", lw=1, alpha=0.6)
axes[1].axhline(np.nanmean(his_diag2), color="teal",   ls="--", lw=1, alpha=0.6)
axes[1].set_title(f"His37 Nδ1 Diagonal Distances  (paper H0: ~7.1 Å)")
axes[1].set_xlabel("Time (ns)")
axes[1].set_ylabel("Distance (Å)")
axes[1].legend()
axes[1].set_xlim(0, times[-1])

fig.suptitle("His37 Tetrad Geometry — AMOEBA Dipole Production")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "his37_tetrad_area.png"), dpi=150)
plt.close()
print("Saved: his37_tetrad_area.png")
print("Done.")
