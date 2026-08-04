"""
AMOEBA NPT production with dipole reporters — H0 (poltype2_equil branch).

Restarts from this folder's own equilibration (traj/equil.rst) rather than a
separate long unrestrained production run, and records total, induced, and
permanent dipoles for every atom at the same interval as the DCD trajectory
(every 10,000 steps = 10 ps).

Output files per run N:
  traj/npt_{N}.dcd                    — trajectory
  traj/npt_{N}.rst                    — restart state
  out/npt_{N}.out                     — energy / temperature / density log
  analysis/total_dipoles_{N}.dat      — total dipole (x y z) per atom per frame
  analysis/induced_dipoles_{N}.dat    — induced dipole per atom per frame
  analysis/permanent_dipoles_{N}.dat  — permanent (lab frame) dipole per atom per frame

Usage:
  python amoeba_npt_dipole.py <counter>
  counter=1  : restart from traj/equil.rst (end of restrained equilibration)
  counter=N  : restart from this folder's traj/npt_{N-1}.rst
"""

import sys
import os
from openmm import *
from openmm.app import *
from openmm.unit import *
from dipolereporter import DipoleReporter, IndDipoleReporter, PermDipoleReporter

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR  = os.path.dirname(SCRIPT_DIR)

STRUCT_PDB  = os.path.join(PARENT_DIR, "amoeba_ready_pdbs", "H0_amoeba_ready.pdb")
LIGAND_XML  = os.path.join(PARENT_DIR, "ligand_amoeba.xml")
POPC_XML    = os.path.join(PARENT_DIR, "popc_poltype2.xml")

# Restart from this folder's own equilibration when counter=1
SEED_RST    = os.path.join(SCRIPT_DIR, "traj", "equil.rst")

for req in [STRUCT_PDB, LIGAND_XML, POPC_XML, SEED_RST]:
    if not os.path.exists(req):
        sys.exit(f"ERROR: required file not found: {req}")

# ── Counter ────────────────────────────────────────────────────────────────────
cnt  = int(sys.argv[1])
pcnt = cnt - 1

if cnt == 1:
    rst = SEED_RST
else:
    rst = os.path.join(SCRIPT_DIR, f"traj/npt_{pcnt}.rst")
    if not os.path.exists(rst):
        sys.exit(f"ERROR: restart file not found: {rst}")

print(f"Run counter : {cnt}")
print(f"Restart from: {rst}")

# ── Simulation parameters ─────────────────────────────────────────────────────
nonbondedCutoff     = 1.2 * nanometers
ewaldErrorTolerance = 0.001
vdwCutoff           = 1.2 * nanometers

dt                  = 0.001 * picoseconds   # 1 fs
temperature         = 310.15 * kelvin
friction            = 1.0 / picosecond
pressure            = 1.0 * atmospheres
barostatInterval    = 25

steps               = 10_000_000            # 10 ns per run
REPORT_INTERVAL     = 10_000               # every 10 ps — matches DCD cadence

platform            = Platform.getPlatformByName("OpenCL")
platformProperties  = {"Precision": "mixed"}

# ── Force field and structure ─────────────────────────────────────────────────
print("Loading force field and structure...")
forcefield = ForceField("amoeba2018.xml", LIGAND_XML, POPC_XML)
pdb = PDBFile(STRUCT_PDB)
print(f"  {pdb.topology.getNumAtoms()} atoms, {pdb.topology.getNumResidues()} residues")

# ── Build system ──────────────────────────────────────────────────────────────
print("Building AMOEBA system...")
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=nonbondedCutoff,
    ewaldErrorTolerance=ewaldErrorTolerance,
    vdwCutoff=vdwCutoff,
    constraints=None,
    rigidWater=False,
    polarization="mutual",
    mutualInducedTargetEpsilon=0.00001,
)
system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))
print(f"  {system.getNumParticles()} particles, {system.getNumForces()} forces")

# ── Integrator ────────────────────────────────────────────────────────────────
integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(1e-6)

simulation = Simulation(pdb.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(pdb.positions)

# ── Load restart state ────────────────────────────────────────────────────────
print(f"Loading restart state...")
with open(rst) as fh:
    simulation.context.setState(XmlSerializer.deserialize(fh.read()))

# ── Reporters ─────────────────────────────────────────────────────────────────
dcd_file  = os.path.join(SCRIPT_DIR, f"traj/npt_{cnt}.dcd")
out_file  = os.path.join(SCRIPT_DIR, f"out/npt_{cnt}.out")
tot_file  = os.path.join(SCRIPT_DIR, f"analysis/total_dipoles_{cnt}.dat")
ind_file  = os.path.join(SCRIPT_DIR, f"analysis/induced_dipoles_{cnt}.dat")
perm_file = os.path.join(SCRIPT_DIR, f"analysis/permanent_dipoles_{cnt}.dat")

simulation.reporters.append(DCDReporter(dcd_file, REPORT_INTERVAL))

simulation.reporters.append(StateDataReporter(
    out_file, REPORT_INTERVAL, totalSteps=steps,
    step=True, progress=True,
    potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
    temperature=True, volume=True, density=True, separator="\t",
))

# AMOEBA-specific dipole reporters
simulation.reporters.append(DipoleReporter(tot_file,  REPORT_INTERVAL))
simulation.reporters.append(IndDipoleReporter(ind_file,  REPORT_INTERVAL))
simulation.reporters.append(PermDipoleReporter(perm_file, REPORT_INTERVAL))

print(f"Dipole reporters writing every {REPORT_INTERVAL} steps ({REPORT_INTERVAL*0.001:.0f} ps)")
print(f"  Total dipoles   -> {tot_file}")
print(f"  Induced dipoles -> {ind_file}")
print(f"  Perm dipoles    -> {perm_file}")

# ── Production ────────────────────────────────────────────────────────────────
print(f"\nRunning {steps:,} steps = {steps*0.001/1000:.0f} ns ...")
simulation.step(steps)

# ── Save restart ──────────────────────────────────────────────────────────────
rst_out = os.path.join(SCRIPT_DIR, f"traj/npt_{cnt}.rst")
state = simulation.context.getState(getPositions=True, getVelocities=True)
with open(rst_out, "w") as fh:
    fh.write(XmlSerializer.serialize(state))
print(f"\nSaved restart: {rst_out}")
print("Done.")
