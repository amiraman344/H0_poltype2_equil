"""
Analysis 10: Channel Water Hydrogen-Bond Network
Water-water and protein-water hydrogen bonds among water confined near the
His37 gate ("channel water"), complementing the RDF/OCF (04/05, His37-water
structure) and ligand H-bond (06, drug-protein) analyses with an explicit
H-bond count/frequency picture of the pore's own hydration network -- the
network thought to support Grotthuss-type proton shuttling through the M2
channel.

Channel water definition (matches script 09's channel/bulk cutoff):
  water O within 8.0 A of the nearest His37 Nd1 atom, evaluated every frame.

H-bond criteria (same as script 06 / the paper):
  donor-acceptor distance <= 4.0 A
  donor-H-acceptor angle  >= 120 deg

Performance design (history, see git/session notes):
  v1: a per-frame dynamic "water within 8 A of His37" selection string was
      far too slow (MDAnalysis rebuilds the selection from scratch every
      frame).
  v2: a whole-system (all 3,469 waters) static-selection H-bond search,
      post-filtered to the channel, was tried next -- fine for the 10 ns
      trajectory but far too slow for 100 ns.
  v3: restricting the whole-100-ns-trajectory search to waters that EVER
      visit the channel over the full 100 ns cut the search space (3,469 ->
      2,120 waters) but still stalled: it ran for 45+ hours of wall clock
      (only ~2 hours of actual CPU time) because of severe memory-swap
      thrashing on this machine, and made no further progress after
      reaching 88% -- killed.
  v4 (this version): process one 10 ns counter (npt_N.dcd) at a time, exactly
      as script 09's induced-dipole analysis does for the same reason. Each
      counter's "ever visits the channel" candidate set is computed and
      searched independently (correctness is unaffected -- a water that
      never enters the channel during a given counter cannot contribute a
      channel H-bond during that counter's frames, regardless of what it
      does in other counters), and only compact per-frame counts / per-
      residue tallies are kept in memory across counters. This bounds peak
      memory to roughly one counter's worth of data instead of the whole
      100 ns trajectory's candidate set at once.

Outputs:
  results/10_channel_hbonds/water_water_hbonds_timeseries.png
  results/10_channel_hbonds/water_water_hbonds.csv
  results/10_channel_hbonds/protein_water_hbonds_frequency.png
  results/10_channel_hbonds/protein_water_hbonds_timeseries.png
  results/10_channel_hbonds/protein_water_hbonds.csv
"""

import os, glob, time
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TOP      = os.path.join(SCRIPT_DIR, "../../amoeba_ready_pdbs/H0_amoeba_ready.pdb")
TRAJ_DIR = os.path.join(SCRIPT_DIR, "../traj")
OUT      = os.path.join(SCRIPT_DIR, "results/10_channel_hbonds")
os.makedirs(OUT, exist_ok=True)

CHANNEL_CUTOFF = 8.0   # A, matches script 09's channel-water definition

def numkey(path):
    return int(os.path.basename(path).split("_")[-1].split(".")[0])

TRAJS = sorted(glob.glob(os.path.join(TRAJ_DIR, "npt_*.dcd")), key=numkey)
if not TRAJS:
    raise FileNotFoundError(f"No npt_*.dcd found in {TRAJ_DIR}")
print(f"Processing {len(TRAJS)} counter(s), one at a time: {[os.path.basename(t) for t in TRAJS]}")

# topology-only universe, just for atom counts/selections shared across counters
u_top = mda.Universe(TOP)
n_atoms_total = len(u_top.atoms)
wat_o_top = u_top.select_atoms("resname HOH and name OH2")
his_top   = u_top.select_atoms("resname HID and name ND1")
n_wat = len(wat_o_top)
print(f"Water O: {n_wat}  His37 ND1: {len(his_top)}")

dt_ps = None
n_frames_per_segment = []
ww_per_frame_chunks = []
n_channel_water_chunks = []
resid_counter = Counter()
pw_per_frame_chunks = []

for seg_i, tf in enumerate(TRAJS, start=1):
    t_seg0 = time.time()
    print(f"\n=== Counter {seg_i}/{len(TRAJS)}: {os.path.basename(tf)} ===")
    useg = mda.Universe(TOP, tf)
    if dt_ps is None:
        dt_ps = useg.trajectory.dt
    n_frames = useg.trajectory.n_frames
    n_frames_per_segment.append(n_frames)

    wat_o   = useg.select_atoms("resname HOH and name OH2")
    his_nd1 = useg.select_atoms("resname HID and name ND1")

    # -- position pass for this counter only --
    dist_to_his = np.empty((n_frames, n_wat), dtype=np.float32)
    n_channel_water = np.empty(n_frames, dtype=np.int32)
    for i, ts in enumerate(useg.trajectory):
        wp = wat_o.positions
        hp = his_nd1.positions
        d = np.linalg.norm(wp[:, None, :] - hp[None, :, :], axis=2).min(axis=1)
        dist_to_his[i] = d
        n_channel_water[i] = np.count_nonzero(d < CHANNEL_CUTOFF)
    useg.trajectory.close()
    n_channel_water_chunks.append(n_channel_water)
    print(f"  mean channel water this counter: {n_channel_water.mean():.1f}")

    # -- candidate restriction: waters ever in-channel during THIS counter --
    ever_channel = dist_to_his.min(axis=0) < CHANNEL_CUTOFF
    n_candidates = int(ever_channel.sum())
    print(f"  candidate waters this counter: {n_candidates} of {n_wat}")

    candidate_residues = wat_o[ever_channel].residues
    candidate_atoms = candidate_residues.atoms
    candidate_index_str = " ".join(str(i) for i in candidate_atoms.indices)
    CANDIDATE_WATER_SEL = f"index {candidate_index_str}"

    atom_to_row = np.full(n_atoms_total, -1, dtype=np.int64)
    atom_to_row[wat_o.indices] = np.arange(n_wat)

    def channel_mask_for(hbonds, water_col):
        frames = hbonds[:, 0].astype(np.int64)
        water_idx = hbonds[:, water_col].astype(np.int64)
        rows = atom_to_row[water_idx]
        valid = rows >= 0
        d = np.full(len(hbonds), np.inf, dtype=np.float32)
        d[valid] = dist_to_his[frames[valid], rows[valid]]
        return valid & (d < CHANNEL_CUTOFF), frames

    # need a fresh universe per HBA call using this counter's file (HBA iterates its own copy)
    useg2 = mda.Universe(TOP, tf)
    hba_ww = HydrogenBondAnalysis(
        universe=useg2,
        donors_sel=f"({CANDIDATE_WATER_SEL}) and name OH2",
        hydrogens_sel=f"({CANDIDATE_WATER_SEL}) and name H1 H2",
        acceptors_sel=f"({CANDIDATE_WATER_SEL}) and name OH2",
        d_a_cutoff=4.0, d_h_a_angle_cutoff=120.0, update_selections=False,
    )
    hba_ww.run(verbose=False)

    ww_this = np.zeros(n_frames)
    hbonds = hba_ww.results.hbonds
    if hbonds is not None and len(hbonds) > 0:
        mask_d, frames = channel_mask_for(hbonds, water_col=1)
        mask_a, _      = channel_mask_for(hbonds, water_col=3)
        channel_mask = mask_d & mask_a
        ww_this = np.bincount(frames[channel_mask], minlength=n_frames).astype(float)
    ww_per_frame_chunks.append(ww_this)
    print(f"  channel water-water H-bonds/frame (mean this counter): {ww_this.mean():.2f}")

    useg3 = mda.Universe(TOP, tf)
    hba_pw1 = HydrogenBondAnalysis(
        universe=useg3, donors_sel="protein", hydrogens_sel="protein and name H*",
        acceptors_sel=f"({CANDIDATE_WATER_SEL}) and name OH2",
        d_a_cutoff=4.0, d_h_a_angle_cutoff=120.0, update_selections=False,
    )
    hba_pw1.run(verbose=False)

    useg4 = mda.Universe(TOP, tf)
    hba_pw2 = HydrogenBondAnalysis(
        universe=useg4, donors_sel=f"({CANDIDATE_WATER_SEL}) and name OH2",
        hydrogens_sel=f"({CANDIDATE_WATER_SEL}) and name H1 H2",
        acceptors_sel="protein and name N* O* S*",
        d_a_cutoff=4.0, d_h_a_angle_cutoff=120.0, update_selections=False,
    )
    hba_pw2.run(verbose=False)

    pw_this = np.zeros(n_frames)

    def accumulate(hba_result, u_ref, water_col, protein_col):
        hbonds = hba_result.results.hbonds
        if hbonds is None or len(hbonds) == 0:
            return
        mask, frames = channel_mask_for(hbonds, water_col=water_col)
        if not mask.any():
            return
        kept_protein_idx = hbonds[mask, protein_col].astype(np.int64)
        resnames = u_ref.atoms.resnames[kept_protein_idx]
        resids   = u_ref.atoms.resids[kept_protein_idx]
        resid_counter.update(f"{rn}{ri}" for rn, ri in zip(resnames, resids))
        global pw_this
        pw_this = pw_this + np.bincount(frames[mask], minlength=n_frames)

    accumulate(hba_pw1, useg3, water_col=3, protein_col=1)
    accumulate(hba_pw2, useg4, water_col=1, protein_col=3)
    pw_per_frame_chunks.append(pw_this)
    print(f"  channel protein-water H-bonds/frame (mean this counter): {pw_this.mean():.2f}")

    del dist_to_his, atom_to_row
    print(f"  counter {seg_i} done in {time.time()-t_seg0:.1f} s")

total_frames = sum(n_frames_per_segment)
times_ns = np.arange(total_frames) * (dt_ps / 1000)
n_channel_water_per_frame = np.concatenate(n_channel_water_chunks)
ww_per_frame = np.concatenate(ww_per_frame_chunks)
pw_per_frame = np.concatenate(pw_per_frame_chunks)

mean_channel_water = n_channel_water_per_frame.mean()
mean_ww_hbonds = ww_per_frame.mean()
hbonds_per_water = mean_ww_hbonds / mean_channel_water if mean_channel_water > 0 else float("nan")
print(f"\n=== TOTAL across {total_frames} frames ===")
print(f"Mean channel water molecules per frame (<{CHANNEL_CUTOFF:.0f} A of His37): {mean_channel_water:.1f}")
print(f"Mean water-water H-bonds per frame (channel-internal): {mean_ww_hbonds:.2f}")
print(f"Mean H-bonds per channel water molecule: {hbonds_per_water:.3f}  "
      f"(bulk water is typically ~3.5 per molecule; a lower value here indicates "
      f"a disrupted/incomplete network near the gate)")

np.savetxt(os.path.join(OUT, "water_water_hbonds.csv"),
           np.column_stack([times_ns, n_channel_water_per_frame, ww_per_frame]),
           header="time_ns,n_channel_water,n_ww_hbonds", delimiter=",", comments="")
print("Saved: water_water_hbonds.csv")

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(times_ns, ww_per_frame, lw=0.5, color="teal", label="Water-water H-bonds (channel)")
ax.axhline(mean_ww_hbonds, color="k", ls="--", lw=1, label=f"Mean: {mean_ww_hbonds:.1f}")
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Number of H-bonds")
ax.set_title(f"Channel-Internal Water-Water H-bonds -- AMOEBA Dipole Production\n"
             f"(waters within {CHANNEL_CUTOFF:.0f} A of nearest His37 Nd1)")
ax.legend()
ax.set_xlim(0, times_ns[-1])
fig.tight_layout()
fig.savefig(os.path.join(OUT, "water_water_hbonds_timeseries.png"), dpi=150)
plt.close()
print("Saved: water_water_hbonds_timeseries.png")

resid_freq = {k: v / total_frames for k, v in resid_counter.items()}
sorted_res = sorted(resid_freq.items(), key=lambda x: -x[1])

with open(os.path.join(OUT, "protein_water_hbonds.csv"), "w") as f:
    f.write("residue,hbond_frequency\n")
    for res, freq in sorted_res:
        f.write(f"{res},{freq:.4f}\n")
print("Saved: protein_water_hbonds.csv")

if sorted_res:
    labels = [r[0] for r in sorted_res[:20]]
    freqs  = [r[1] for r in sorted_res[:20]]
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(labels, freqs, color="mediumseagreen", edgecolor="k", lw=0.5)
    ax.set_xlabel("Protein Residue")
    ax.set_ylabel("H-bond Frequency (fraction of frames)")
    ax.set_title(f"Protein - Channel-Water H-bond Frequency\n"
                 f"(waters within {CHANNEL_CUTOFF:.0f} A of nearest His37 Nd1)")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "protein_water_hbonds_frequency.png"), dpi=150)
    plt.close()
    print("Saved: protein_water_hbonds_frequency.png")
    print("Top 10 protein-channel-water H-bond partners:")
    for r, fr in sorted_res[:10]:
        print(f"  {r}: {fr*100:.1f}% of frames")
else:
    print("No protein-channel-water H-bonds detected.")

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(times_ns, pw_per_frame, lw=0.5, color="mediumseagreen")
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Number of H-bonds")
ax.set_title("Protein-Channel-Water H-bonds per Frame -- AMOEBA Dipole Production")
ax.set_xlim(0, times_ns[-1])
fig.tight_layout()
fig.savefig(os.path.join(OUT, "protein_water_hbonds_timeseries.png"), dpi=150)
plt.close()
print("Saved: protein_water_hbonds_timeseries.png")
print("Done.")
