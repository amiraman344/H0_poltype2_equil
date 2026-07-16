"""
Unrestrained NPT production — 10 ns.

Starts from traj/equil.rst (output of equil_npt.py).
No position restraints. All AMOEBA forces active.

Output:
  traj/npt_1.dcd / npt_1.rst  (counter=1 → 10 ns)
  traj/npt_N.dcd / npt_N.rst  (counter=N → restart from npt_{N-1}.rst)
"""

import sys, os
from openmm import *
from openmm.app import *
from openmm.unit import *

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
COMPLEX_PDB = os.path.join(SCRIPT_DIR, "complex.pdb")
LIGAND_XML  = os.path.join(SCRIPT_DIR, "ligand_amoeba.xml")
POPC_XML    = os.path.join(SCRIPT_DIR, "popc_poltype2.xml")

for req in [COMPLEX_PDB, LIGAND_XML, POPC_XML]:
    if not os.path.exists(req):
        sys.exit(f"ERROR: {req} not found.")

cnt  = int(sys.argv[1])
pcnt = cnt - 1

# counter=1 restarts from equilibration output; counter>1 from previous production run
if cnt == 1:
    rst = os.path.join(SCRIPT_DIR, "traj/equil.rst")
else:
    rst = os.path.join(SCRIPT_DIR, f"traj/npt_{pcnt}.rst")

if not os.path.exists(rst):
    sys.exit(f"ERROR: restart file not found: {rst}")

for d in ["traj", "out"]:
    os.makedirs(os.path.join(SCRIPT_DIR, d), exist_ok=True)

# ── Simulation parameters ──────────────────────────────────────────────────────
nonbondedCutoff     = 1.2 * nanometers
ewaldErrorTolerance = 0.001
vdwCutoff           = 1.2 * nanometers

dt                  = 0.001 * picoseconds
temperature         = 310.15 * kelvin
friction            = 1.0 / picosecond
pressure            = 1.0 * atmospheres
barostatInterval    = 25

steps               = 10_000_000           # 10 ns per run
reportInterval      = 10_000              # save frame every 10 ps

platform = Platform.getPlatformByName("OpenCL")
platformProperties = {"Precision": "mixed"}

# ── Load system ────────────────────────────────────────────────────────────────
print("Loading force field and structure...")
forcefield = ForceField("amoeba2018.xml", LIGAND_XML, POPC_XML)
pdb = PDBFile(COMPLEX_PDB)
print(f"  {pdb.topology.getNumAtoms()} atoms, {pdb.topology.getNumResidues()} residues")

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

# ── Integrator and simulation ──────────────────────────────────────────────────
integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(1e-6)

simulation = Simulation(pdb.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(pdb.positions)

# ── Load restart ───────────────────────────────────────────────────────────────
print(f"Restarting from {rst}")
with open(rst) as fh:
    simulation.context.setState(XmlSerializer.deserialize(fh.read()))

# ── Reporters ──────────────────────────────────────────────────────────────────
simulation.reporters.append(DCDReporter(
    os.path.join(SCRIPT_DIR, f"traj/npt_{cnt}.dcd"), reportInterval
))
simulation.reporters.append(StateDataReporter(
    os.path.join(SCRIPT_DIR, f"out/npt_{cnt}.out"),
    reportInterval, totalSteps=steps,
    step=True, progress=True,
    potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
    temperature=True, volume=True, density=True, separator="\t",
))

# ── Production MD ──────────────────────────────────────────────────────────────
print(f"Running production MD ({steps} steps = {steps*0.001:.0f} ps)...")
simulation.step(steps)

state = simulation.context.getState(getPositions=True, getVelocities=True)
rst_out = os.path.join(SCRIPT_DIR, f"traj/npt_{cnt}.rst")
with open(rst_out, "w") as fh:
    fh.write(XmlSerializer.serialize(state))
print(f"Saved restart: {rst_out}")
