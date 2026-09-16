"""
Analysis 09: AMOEBA Induced Dipole Analysis
UNIQUE TO AMOEBA — not possible with CHARMM/TIP3P/SWM4.

Reads the pre-saved dipole dat files from ../analysis/:
  induced_dipoles_{N}.dat   — induced dipole (x y z) per atom per frame

File format: each line is "x y z" for one atom in atomic units (e·bohr).
Frames are written sequentially: N_atoms lines per frame.
Conversion: 1 e·nm = 48.0321 Debye.

Measures:
  1. Water induced dipole magnitude: channel water (<8 Å from His37) vs bulk
  2. Water induced dipole vs distance from His37 Nδ1
  3. Per-residue protein induced dipole magnitudes
  4. Ligand induced dipole over time
  5. Polarization enhancement inside channel vs bulk

Outputs:
  results/induced_dipole_water_vs_distance.png
  results/induced_dipole_channel_vs_bulk.png
  results/induced_dipole_protein.png
  results/induced_dipole_ligand.png
  results/induced_dipoles_water.csv
  results/induced_dipoles_protein.csv
  results/induced_dipoles_ligand.csv

NOTE (100 ns version): rewritten to process one npt_N.dcd / induced_dipoles_N.dat
counter at a time instead of concatenating all counters into memory at once.
With 10 counters x ~19,723 atoms x 1000 frames, loading everything simultaneously
(as the original 10 ns script did) would require several GB of RAM per array and
risks OOM on shared/interactive nodes. Segment-at-a-time processing bounds peak
memory to roughly one counter's worth of data while producing numerically
identical results. Also fixes a trajectory/dipole-file ordering bug: plain
string sorted() on "npt_10.dcd" et al. sorts lexicographically (1, 10, 2, 3, ...)
rather than numerically, which would silently scramble frame order once a
double-digit counter exists.
"""

import os, glob
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
TOP         = os.path.join(SCRIPT_DIR, "../../amoeba_ready_pdbs/H0_amoeba_ready.pdb")
TRAJ_DIR    = os.path.join(SCRIPT_DIR, "../traj")
DIPOLE_DIR  = os.path.join(SCRIPT_DIR, "../analysis")
OUT         = os.path.join(SCRIPT_DIR, "results/09_induced_dipoles")
os.makedirs(OUT, exist_ok=True)

# Conversion: OpenMM DipoleReporter stores raw values in e·nm = 48.0321 D
EBOHR_TO_DEBYE = 48.0321

CHANNEL_CUTOFF_A = 8.0
BULK_CUTOFF_A    = 15.0


def numkey(path):
    """Numeric sort key: 'npt_10.dcd' -> 10, 'induced_dipoles_2.dat' -> 2."""
    return int(os.path.basename(path).split("_")[-1].split(".")[0])


traj_files = sorted(glob.glob(os.path.join(TRAJ_DIR, "npt_*.dcd")), key=numkey)
dat_files  = sorted(glob.glob(os.path.join(DIPOLE_DIR, "induced_dipoles_*.dat")), key=numkey)

if not traj_files:
    raise FileNotFoundError(f"No npt_*.dcd files found in {TRAJ_DIR}")
if not dat_files:
    raise FileNotFoundError(
        f"No induced_dipoles_*.dat files found in {DIPOLE_DIR}\n"
        f"Run the dipole production simulation first."
    )

traj_counters = [numkey(p) for p in traj_files]
dat_counters  = [numkey(p) for p in dat_files]
if traj_counters != dat_counters:
    raise RuntimeError(
        f"Trajectory counters {traj_counters} do not match dipole-file counters "
        f"{dat_counters} — refusing to guess an alignment."
    )
print(f"Segments to process (in order): {traj_counters}")

# ── Topology-only universe for atom selections ────────────────────────────────
u_top = mda.Universe(TOP)
N_ATOMS = len(u_top.atoms)
print(f"System: {N_ATOMS} atoms")

wat_o    = u_top.select_atoms("resname HOH and name OH2")
his_nd1  = u_top.select_atoms("resname HID and name ND1")
ligand   = u_top.select_atoms("resname 308")
protein  = u_top.select_atoms("protein")

wat_idx  = wat_o.indices
his_idx  = his_nd1.indices
lig_idx  = ligand.indices

print(f"Water O: {len(wat_idx)}  His ND1: {len(his_idx)}  Ligand: {len(lig_idx)}")


def load_dipole_dat(dat_path, n_atoms):
    """Load one dipole dat file -> (n_frames, n_atoms, 3) array, in Debye, float32."""
    raw = np.loadtxt(dat_path, dtype=np.float32)          # (n_frames*n_atoms, 3)
    n_frames = raw.shape[0] // n_atoms
    dipoles = raw[: n_frames * n_atoms].reshape(n_frames, n_atoms, 3)
    return dipoles * np.float32(EBOHR_TO_DEBYE)


# ── Running accumulators across segments ──────────────────────────────────────
dt_ps = None
water_dist_chunks = []   # per-frame (N_wat,) distance arrays
water_mu_chunks   = []   # per-frame (N_wat,) induced-mag arrays
lig_mean_chunks   = []   # per-segment (n_use,) arrays
lig_max_chunks    = []
resid_stat = {}          # resid -> [sum, count, resname]
n_use_total = 0

for seg_i, (tf, df) in enumerate(zip(traj_files, dat_files), start=1):
    print(f"\n--- Segment {seg_i}/{len(traj_files)}: "
          f"{os.path.basename(tf)} + {os.path.basename(df)} ---")

    useg = mda.Universe(TOP, tf)
    if dt_ps is None:
        dt_ps = useg.trajectory.dt
    n_traj_frames = useg.trajectory.n_frames

    # Selections must be made on THIS universe (useg), not the topology-only
    # u_top, so that .positions tracks useg's currently-loaded trajectory frame.
    wat_o_seg   = useg.select_atoms("resname HOH and name OH2")
    his_nd1_seg = useg.select_atoms("resname HID and name ND1")

    wat_pos = np.empty((n_traj_frames, len(wat_idx), 3), dtype=np.float32)
    his_pos = np.empty((n_traj_frames, len(his_idx), 3), dtype=np.float32)
    for fi, ts in enumerate(useg.trajectory):
        wat_pos[fi] = wat_o_seg.positions
        his_pos[fi] = his_nd1_seg.positions
    useg.trajectory.close()

    print(f"  Loading {os.path.basename(df)}...")
    dip = load_dipole_dat(df, N_ATOMS)      # (n_dip_frames, N_ATOMS, 3)
    n_dip_frames = dip.shape[0]

    n_use = min(n_traj_frames, n_dip_frames)
    if n_traj_frames != n_dip_frames:
        print(f"  WARNING: frame count mismatch — traj={n_traj_frames} "
              f"dipole={n_dip_frames}; using first {n_use} frames of each")

    induced_mag = np.linalg.norm(dip[:n_use], axis=2)   # (n_use, N_ATOMS)
    del dip

    for fi in range(n_use):
        wp = wat_pos[fi]
        hp = his_pos[fi]
        if len(hp) > 0:
            d = np.min(
                np.linalg.norm(wp[:, None, :] - hp[None, :, :], axis=2), axis=1
            )
        else:
            d = np.full(len(wat_idx), np.nan, dtype=np.float32)
        water_dist_chunks.append(d)
        water_mu_chunks.append(induced_mag[fi, wat_idx])

        if fi % 200 == 0:
            print(f"    frame {fi}/{n_use}")

    lig_mean_chunks.append(induced_mag[:, lig_idx].mean(axis=1))
    lig_max_chunks.append(induced_mag[:, lig_idx].max(axis=1))

    for res in protein.residues:
        idx = res.atoms.indices
        s = float(induced_mag[:, idx].sum())
        c = induced_mag[:, idx].size
        if res.resid not in resid_stat:
            resid_stat[res.resid] = [0.0, 0, res.resname]
        resid_stat[res.resid][0] += s
        resid_stat[res.resid][1] += c

    n_use_total += n_use
    print(f"  frames used this segment: {n_use}  (running total: {n_use_total})")

    del induced_mag, wat_pos, his_pos

# ── Assemble frame-major arrays (same layout as the original single-pass script) ──
dist_col = np.concatenate(water_dist_chunks)     # (n_use_total * N_wat,)
mu_col   = np.concatenate(water_mu_chunks)
frame_col = np.repeat(np.arange(n_use_total), len(wat_idx))
records_water = np.column_stack([frame_col, dist_col, mu_col])

lig_mean = np.concatenate(lig_mean_chunks)       # (n_use_total,)
lig_max  = np.concatenate(lig_max_chunks)

prot_records = [(resid, resname, s / c) for resid, (s, c, resname) in resid_stat.items()]

print(f"\nTotal frames analyzed: {n_use_total}")

# ── Save CSVs ─────────────────────────────────────────────────────────────────
np.savetxt(os.path.join(OUT, "induced_dipoles_water.csv"),
           records_water,
           header="frame,dist_to_HisND1_A,induced_D",
           delimiter=",", comments="")

with open(os.path.join(OUT, "induced_dipoles_protein.csv"), "w") as f:
    f.write("resid,resname,mean_induced_D\n")
    for resid, resname, mu in sorted(prot_records, key=lambda x: x[0]):
        f.write(f"{resid},{resname},{mu:.5f}\n")

times_ns = np.arange(n_use_total) * (dt_ps / 1000)
np.savetxt(os.path.join(OUT, "induced_dipoles_ligand.csv"),
           np.column_stack([times_ns, lig_mean, lig_max]),
           header="time_ns,mean_induced_D,max_induced_D",
           delimiter=",", comments="")

print("Saved: induced_dipoles_water.csv, induced_dipoles_protein.csv, induced_dipoles_ligand.csv")

# ── Bulk water baseline (>15 Å from His37) ───────────────────────────────────
channel_mu = mu_col[dist_col < CHANNEL_CUTOFF_A]
bulk_mu    = mu_col[dist_col > BULK_CUTOFF_A]
print(f"\nChannel water (<{CHANNEL_CUTOFF_A:.0f} Å from His37) mean |mu_ind| = {channel_mu.mean():.4f} D")
print(f"Bulk water    (>{BULK_CUTOFF_A:.0f} Å from His37) mean |mu_ind| = {bulk_mu.mean():.4f} D")
print(f"Ratio (channel/bulk): {channel_mu.mean()/bulk_mu.mean():.3f}x")

# ── Plot 1: water dipole vs distance from His37 ───────────────────────────────
dist_bins   = np.linspace(0, 20, 40)
bin_centers = 0.5 * (dist_bins[:-1] + dist_bins[1:])
bin_idx     = np.digitize(dist_col, dist_bins)
bin_mean    = np.array([mu_col[bin_idx == k].mean() if (bin_idx == k).any() else np.nan
                        for k in range(1, len(dist_bins))])
bin_std     = np.array([mu_col[bin_idx == k].std()  if (bin_idx == k).any() else 0
                        for k in range(1, len(dist_bins))])

fig, ax = plt.subplots(figsize=(8, 5))
ax.fill_between(bin_centers, bin_mean - bin_std, bin_mean + bin_std,
                alpha=0.3, color="steelblue")
ax.plot(bin_centers, bin_mean, color="steelblue", lw=1.5, label="Mean |mu_ind|")
bulk_base = np.nanmean(bin_mean[bin_centers > 15])
ax.axhline(bulk_base, color="gray", ls="--", lw=1,
           label=f"Bulk baseline: {bulk_base:.4f} D")
ax.set_xlabel("Distance to nearest His37 Nδ1 (Å)")
ax.set_ylabel("|mu_ind| water O (Debye)")
ax.set_title("AMOEBA Water Induced Dipole vs Distance from His37\n"
             "(Novel AMOEBA analysis — not in Gődény et al. 2025)")
ax.legend()
ax.set_xlim(0, 20)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "induced_dipole_water_vs_distance.png"), dpi=150)
plt.close()
print("Saved: induced_dipole_water_vs_distance.png")

# ── Plot 2: channel vs bulk histogram ────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5))
ax.hist(channel_mu, bins=50, alpha=0.6, color="crimson", density=True,
        label=f"Channel water (<8 Å from His)  n={len(channel_mu):,}")
ax.hist(bulk_mu,    bins=50, alpha=0.6, color="steelblue", density=True,
        label=f"Bulk water (>15 Å from His)  n={len(bulk_mu):,}")
ax.set_xlabel("|mu_ind| water O (Debye)")
ax.set_ylabel("Probability density")
ax.set_title("AMOEBA Induced Dipole: Channel vs Bulk Water")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "induced_dipole_channel_vs_bulk.png"), dpi=150)
plt.close()
print("Saved: induced_dipole_channel_vs_bulk.png")

# ── Plot 3: per-residue protein induced dipoles ───────────────────────────────
prot_sorted = sorted(prot_records, key=lambda x: x[0])
resids   = [r[0] for r in prot_sorted]
res_mu   = [r[2] for r in prot_sorted]
resnames = [r[1] for r in prot_sorted]
colors   = ["crimson" if rn == "HID" else "steelblue" for rn in resnames]

fig, ax = plt.subplots(figsize=(10, 4))
ax.bar(resids, res_mu, color=colors, edgecolor="none", width=0.8)
his_patch  = plt.Rectangle((0, 0), 1, 1, color="crimson", label="His37 (HID)")
prot_patch = plt.Rectangle((0, 0), 1, 1, color="steelblue", label="Other protein")
ax.legend(handles=[his_patch, prot_patch])
ax.set_xlabel("Residue ID")
ax.set_ylabel("Mean |mu_ind| (Debye)")
ax.set_title("Per-Residue Protein Induced Dipole — AMOEBA\n"
             "(Novel analysis: polarization at each residue)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "induced_dipole_protein.png"), dpi=150)
plt.close()
print("Saved: induced_dipole_protein.png")

# ── Plot 4: ligand induced dipole over time ───────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(times_ns, lig_mean, color="darkorange", lw=1.2, label="Mean |mu_ind|")
ax.fill_between(times_ns, lig_mean * 0.9, lig_max, alpha=0.3, color="darkorange")
ax.set_xlabel("Time (ns)")
ax.set_ylabel("|mu_ind| (Debye)")
ax.set_title("Ligand (Amantadine) Induced Dipole — AMOEBA Dipole Production")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(OUT, "induced_dipole_ligand.png"), dpi=150)
plt.close()
print("Saved: induced_dipole_ligand.png")
print("\nDone.")
