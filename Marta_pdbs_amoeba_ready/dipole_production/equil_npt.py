"""
Restrained NPT equilibration — 5 ns — starting point for dipole production.
Marta's H0 system (apo M2 channel, no bound ligand).

Harmonic position restraints are applied to:
  - Protein Cα atoms          k = 1000 kJ/mol/nm²
  - Protein backbone (N,C,O)  k =  500 kJ/mol/nm²
  - POPC phosphorus (P)       k =  500 kJ/mol/nm²

Restraints hold the protein and membrane headgroups in place while
lipid tails and water relax around them at 310.15 K, 1 atm.

Starting structure: ../H0_amoeba_ready.pdb (converted from Marta's raw H0.pdb
by ../prepare_marta_amoeba.py — CRYST1 box, renamed residues, CONECT records,
and each chain's C-terminal Leu46 fused with an N-methylamide cap).

Differs from the H0_poltype2_equil/dipole_production reference this is based
on in two ways:
  1. No ligand — this system is apo, so no ligand force-field XML is loaded.
  2. Loads ../leu_cterm_amide_cap.xml alongside amoeba2018.xml — the fused
     Leu46+cap residue has no template in vanilla amoeba2018.xml.

Output:
  traj/equil.dcd   — trajectory (1 frame / 50 ps)
  traj/equil.rst   — restart, seeds amoeba_npt_dipole.py counter=1
"""

import sys, os
import numpy as np
from openmm import *
from openmm.app import *
from openmm.unit import *

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR  = os.path.dirname(SCRIPT_DIR)

STRUCT_PDB  = os.path.join(PARENT_DIR, "H0_amoeba_ready.pdb")
LEU_CAP_XML = os.path.join(PARENT_DIR, "leu_cterm_amide_cap.xml")
POPC_XML    = "/tab3/aman/AMEOBA_project/poltype2_popc/H0_poltype2_equil/popc_poltype2.xml"

for req in [STRUCT_PDB, LEU_CAP_XML, POPC_XML]:
    if not os.path.exists(req):
        sys.exit(f"ERROR: {req} not found.")

for d in ["traj", "out"]:
    os.makedirs(os.path.join(SCRIPT_DIR, d), exist_ok=True)

# ── Simulation parameters ──────────────────────────────────────────────────────
nonbondedCutoff     = 1.2 * nanometers
ewaldErrorTolerance = 0.001
vdwCutoff           = 1.2 * nanometers

dt                  = 0.001 * picoseconds   # 1 fs
temperature         = 310.15 * kelvin
friction            = 1.0 / picosecond
pressure            = 1.0 * atmospheres
barostatInterval    = 25

equilSteps          = 5_000_000             # 5 ns
miniSteps           = 1_000                 # energy minimisation iterations
reportInterval      = 50_000               # save frame every 50 ps

# Restraint spring constants
K_CA       = 1000 * kilojoules_per_mole / nanometer**2   # protein Cα
K_BB       =  500 * kilojoules_per_mole / nanometer**2   # backbone N, C, O
K_P        =  500 * kilojoules_per_mole / nanometer**2   # POPC phosphorus

platform = Platform.getPlatformByName("OpenCL")
platformProperties = {"Precision": "mixed"}

# ── Load system ────────────────────────────────────────────────────────────────
print("Loading force field and structure...")
forcefield = ForceField("amoeba2018.xml", LEU_CAP_XML, POPC_XML)
pdb = PDBFile(STRUCT_PDB)
print(f"  {pdb.topology.getNumAtoms()} atoms, {pdb.topology.getNumResidues()} residues")
print(f"  Box: {pdb.topology.getPeriodicBoxVectors()}")

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

# ── Position restraints ────────────────────────────────────────────────────────
print("Adding position restraints...")
restraint = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
restraint.addGlobalParameter("k", 1.0)   # placeholder, overridden per-particle via expression
restraint.addPerParticleParameter("k")
restraint.addPerParticleParameter("x0")
restraint.addPerParticleParameter("y0")
restraint.addPerParticleParameter("z0")

# Build lookup: atom index -> (x0, y0, z0) in nm
positions_nm = pdb.positions.in_units_of(nanometers)

n_ca = n_bb = n_p = 0
for atom in pdb.topology.atoms():
    res  = atom.residue
    name = atom.name
    pos  = positions_nm[atom.index]
    x0, y0, z0 = pos.x, pos.y, pos.z

    if atom.element is not None and res.name not in ("HOH", "POPC"):
        # protein atoms (no ligand in this system)
        if name == "CA":
            restraint.addParticle(atom.index, [K_CA.value_in_unit(kilojoules_per_mole/nanometer**2), x0, y0, z0])
            n_ca += 1
        elif name in ("N", "C", "O"):
            restraint.addParticle(atom.index, [K_BB.value_in_unit(kilojoules_per_mole/nanometer**2), x0, y0, z0])
            n_bb += 1
    elif res.name == "POPC" and name == "P":
        restraint.addParticle(atom.index, [K_P.value_in_unit(kilojoules_per_mole/nanometer**2), x0, y0, z0])
        n_p += 1

system.addForce(restraint)
print(f"  Restrained: {n_ca} Cα (k={K_CA}), {n_bb} backbone N/C/O (k={K_BB}), {n_p} POPC-P (k={K_P})")
print(f"  {system.getNumParticles()} particles, {system.getNumForces()} forces")

# ── Integrator and simulation ──────────────────────────────────────────────────
integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(1e-6)

simulation = Simulation(pdb.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(pdb.positions)

# ── Minimise ───────────────────────────────────────────────────────────────────
print(f"Minimising energy ({miniSteps} iterations)...")
simulation.minimizeEnergy(maxIterations=miniSteps)

# ── Reporters ──────────────────────────────────────────────────────────────────
simulation.reporters.append(DCDReporter(
    os.path.join(SCRIPT_DIR, "traj/equil.dcd"), reportInterval
))
simulation.reporters.append(StateDataReporter(
    os.path.join(SCRIPT_DIR, "out/equil.out"),
    reportInterval, totalSteps=equilSteps,
    step=True, progress=True,
    potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
    temperature=True, volume=True, density=True, separator="\t",
))

# ── Equilibration ──────────────────────────────────────────────────────────────
simulation.context.setVelocitiesToTemperature(temperature)
print(f"Running restrained NPT equilibration ({equilSteps} steps = {equilSteps*0.001:.0f} ps)...")
simulation.step(equilSteps)

state = simulation.context.getState(getPositions=True, getVelocities=True)
rst_out = os.path.join(SCRIPT_DIR, "traj/equil.rst")
with open(rst_out, "w") as fh:
    fh.write(XmlSerializer.serialize(state))
print(f"Saved restart: {rst_out}")
