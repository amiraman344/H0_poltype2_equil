"""
Convert Marta's raw H0-H4/H2_trans PDBs to AMOEBA-compatible format.

Fixes applied (parallel to check_nutralized_system/prepare_all_amoeba.py
used for our own H0-H4 systems):
  1. Residue renames: TIP3 -> HOH, HSD -> HID, HSP -> HIS
  2. C-terminal cap fix: each chain's Leu46 ends in a bare primary amide
     (-C(=O)NH2, atoms NT/HT1/HT2) which has no AMOEBA2018 template.
     Extended in place (same residue, resName "LEU", resSeq 46 -- NOT a
     separate NME residue) into a full N-methylamide cap:
       - NT  -> kept as-is (position unchanged)
       - HT1 -> kept as-is (position unchanged)
       - HT2 -> discarded; direction used to place the new CT (methyl C)
       - CT + HT2/HT3/HT4 (new methyl group) built via tetrahedral geometry
     Requires the custom residue template leu_cterm_amide_cap.xml (in this
     same folder), loaded alongside amoeba2018.xml, which declares this
     25-atom fused Leu+cap residue reusing amoeba2018's own atom types
     (LEU's for atoms 0-18, NME's for the cap atoms) -- no bare "LEU"
     name/atom-count match exists in vanilla amoeba2018.xml for this.
  3. CRYST1 box record (from coordinate extents + 2 A padding)
  4. CONECT records for POPC (133 bonds/molecule), from popc_poltype2.xml
     bond topology -- no ligand present in these systems, so no ligand
     CONECT needed.
  5. All atom serials renumbered sequentially; TER chain-break records
     are preserved in place (dropping them merges neighboring chains'
     bond topology and breaks System building).

Output validated per-file by building an OpenMM System and running a
short energy minimization (same bar as the original amoeba_ready_pdbs).
"""

import os
import sys
import math
import xml.etree.ElementTree as ET
from collections import defaultdict

IN_DIR  = "/tab3/aman/AMEOBA_project/poltype2_popc/Marta_pdbs"
OUT_DIR = "/tab3/aman/AMEOBA_project/poltype2_popc/Marta_pdbs_amoeba_ready"
POPC_XML = "/tab3/aman/AMEOBA_project/poltype2_popc/H0_poltype2_equil/popc_poltype2.xml"

RENAME = {"TIP3": "HOH", "HSD": "HID", "HSP": "HIS"}

FILES = {
    "H0": "H0.pdb",
    "H1": "H1.pdb",
    "H2": "H2.pdb",
    "H2_trans": "H2_trans.pdb",
    "H3": "H3.pdb",
    "H4": "H4.pdb",
}

CHAINS = ["PROA", "PROB", "PROC", "PROD"]

# ── PDB fixed-column helpers ───────────────────────────────────────────────
def parse_atom(line):
    return dict(
        serial=int(line[6:11]),
        name=line[12:16].strip(),
        altloc=line[16:17],
        resname=line[17:21].strip(),
        chain=line[21:22],
        resnum=int(line[22:26]),
        icode=line[26:27],
        x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
        occ=float(line[54:60]) if line[54:60].strip() else 1.0,
        beta=float(line[60:66]) if line[60:66].strip() else 0.0,
        segid=line[72:76].strip(),
        element=line[76:78].strip(),
    )

def fmt_name(name):
    return name[:4] if len(name) >= 4 else (" " + name).ljust(4)

def fmt_atom(a):
    return (f"ATOM  {a['serial']:5d} {fmt_name(a['name'])}{a['altloc']}"
            f"{a['resname'].ljust(4)[:4]}{a['chain']}{a['resnum']:4d}{a['icode']}   "
            f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}{a['occ']:6.2f}{a['beta']:6.2f}"
            f"      {a['segid'].ljust(4)[:4]}{a['element'].rjust(2)}\n")

# ── Methyl-group geometry ───────────────────────────────────────────────────
def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def add(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def scale(a, s): return (a[0]*s, a[1]*s, a[2]*s)
def norm(a): return math.sqrt(a[0]**2+a[1]**2+a[2]**2)
def normalize(a):
    n = norm(a)
    return (a[0]/n, a[1]/n, a[2]/n)
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def dot(a, b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def build_methyl(p_n, p_ht2):
    """Given N position and the old (discarded) second-H position, build a
    CH3 carbon + 3 hydrogens with proper tetrahedral geometry."""
    axis = normalize(sub(p_ht2, p_n))          # N -> C direction
    p_ch3 = add(p_n, scale(axis, 1.47))        # N-C bond length

    arbitrary = (1.0, 0.0, 0.0) if abs(axis[0]) < 0.9 else (0.0, 1.0, 0.0)
    u = normalize(cross(axis, arbitrary))
    v = cross(axis, u)

    hs = []
    for i in range(3):
        phi = math.radians(i * 120.0)
        dir_h = add(scale(axis, 1.0/3.0),
                    scale(add(scale(u, math.cos(phi)), scale(v, math.sin(phi))),
                          math.sqrt(8.0)/3.0))
        dir_h = normalize(dir_h)
        hs.append(add(p_ch3, scale(dir_h, 1.09)))
    return p_ch3, hs

# ── POPC bond topology (for CONECT) ─────────────────────────────────────────
def parse_xml_bonds(xml_path, resname):
    tree = ET.parse(xml_path)
    for res in tree.findall('.//Residue'):
        if res.get('name') == resname:
            atoms = [a.get('name') for a in res.findall('Atom')]
            bonds = [(int(b.get('from')), int(b.get('to')))
                     for b in res.findall('Bond')]
            return atoms, bonds
    raise RuntimeError(f"Residue '{resname}' not found in {xml_path}")

def conect_from_serial_pairs(serial_pairs):
    bonds_dict = defaultdict(list)
    for s1, s2 in serial_pairs:
        bonds_dict[s1].append(s2)
        bonds_dict[s2].append(s1)
    lines = []
    for s1 in sorted(bonds_dict):
        partners = sorted(bonds_dict[s1])
        for start in range(0, len(partners), 4):
            chunk = partners[start:start+4]
            lines.append("CONECT" + f"{s1:5d}" +
                         "".join(f"{p:5d}" for p in chunk) + "\n")
    return lines

def make_conect(serial_by_name, xml_atoms, xml_bonds):
    pairs = []
    for i, j in xml_bonds:
        a1, a2 = xml_atoms[i], xml_atoms[j]
        if a1 not in serial_by_name or a2 not in serial_by_name:
            continue
        pairs.append((serial_by_name[a1], serial_by_name[a2]))
    return conect_from_serial_pairs(pairs)

popc_atoms, popc_bonds = parse_xml_bonds(POPC_XML, "POPC")
print(f"POPC template: {len(popc_atoms)} atoms, {len(popc_bonds)} bonds\n")

os.makedirs(OUT_DIR, exist_ok=True)

# ── Process each system ─────────────────────────────────────────────────────
target = sys.argv[1] if len(sys.argv) > 1 else None

for sysname, fname in FILES.items():
    if target and sysname != target:
        continue

    in_path = os.path.join(IN_DIR, fname)
    out_path = os.path.join(OUT_DIR, f"{sysname}_amoeba_ready.pdb")
    print(f"[{sysname}] {fname} -> {os.path.basename(out_path)}")

    with open(in_path) as f:
        all_lines = f.readlines()

    # Apply residue renames (in-place on ATOM/HETATM lines) and replace
    # each chain's cap (NT/HT1/HT2 on res 46) with a full amide-cap in place.
    # TER lines (chain-break markers) are passed through untouched, in
    # their original sequential position, so chain topology is preserved.
    tokens = []   # list of atom-dicts and "TER" strings, in final output order
    renamed_count = defaultdict(int)
    n_caps_fixed = 0
    cap_junctions = []  # (c_atom, nam, ham, cam, ham1, ham2, ham3) dicts, for CONECT

    i = 0
    n = len(all_lines)
    while i < n:
        line = all_lines[i]
        if line.startswith("TER"):
            tokens.append("TER")
            i += 1
            continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            i += 1
            continue

        a = parse_atom(line)
        if a["resname"] in RENAME:
            renamed_count[a["resname"]] += 1
            a["resname"] = RENAME[a["resname"]]

        if a["segid"] in CHAINS and a["resnum"] == 46 and a["name"] == "NT":
            b = parse_atom(all_lines[i+1]); assert b["name"] == "HT1"
            c = parse_atom(all_lines[i+2]); assert c["name"] == "HT2"

            p_nt  = (a["x"], a["y"], a["z"])
            p_ht1 = (b["x"], b["y"], b["z"])
            p_old_ht2 = (c["x"], c["y"], c["z"])
            p_ct, hs = build_methyl(p_nt, p_old_ht2)

            template = dict(altloc=" ", chain=a["chain"], icode=" ",
                             occ=1.0, beta=0.0, segid=a["segid"], resname="LEU", resnum=46)

            def mk(name, pos, element):
                d = dict(template)
                d.update(name=name, x=pos[0], y=pos[1], z=pos[2], element=element, serial=0)
                return d

            # Same residue (LEU 46) as the atoms already emitted for this
            # chain -- these just extend it with the amide-cap atoms.
            #
            # NOTE: names deliberately avoid "HT1/HT2/HT3" (and "N"/"H"
            # collisions with the backbone atoms already in this residue).
            # OpenMM's pdbNames.xml aliases HT1/HT2/HT3 -> H/H2/H3 *for any
            # protein residue*, and residues.xml then bonds H2/H3 directly
            # to the residue's own "N" (the standard Amber N-terminal NH3+
            # pattern) -- silently wiring these cap hydrogens to the wrong
            # nitrogen (the backbone N, not the cap's own N) instead of
            # leaving them for distance-based bonding.
            # OpenMM's standard-bonds builder only wires up atoms it
            # recognizes by name within a known residue (LEU); it silently
            # leaves these new cap atoms bond-less otherwise, so their
            # bonds are supplied explicitly via CONECT below. Order in
            # file so far for this residue: ...,C,O -- tokens[-2] is "C".
            c_atom = tokens[-2]
            assert c_atom["name"] == "C", c_atom

            nam  = mk("NAM",  p_nt,  "N")
            ham  = mk("HAM",  p_ht1, "H")
            cam  = mk("CAM",  p_ct,  "C")
            ham1 = mk("HAM1", hs[0], "H")
            ham2 = mk("HAM2", hs[1], "H")
            ham3 = mk("HAM3", hs[2], "H")
            tokens.append(nam); tokens.append(ham); tokens.append(cam)
            tokens.append(ham1); tokens.append(ham2); tokens.append(ham3)

            cap_junctions.append((c_atom, nam, ham, cam, ham1, ham2, ham3))

            n_caps_fixed += 1
            i += 3
            continue

        tokens.append(a)
        i += 1

    assert n_caps_fixed == 4, f"expected 4 chain caps fixed, got {n_caps_fixed}"

    # Renumber serials sequentially (TER tokens are not numbered)
    new_atoms = [t for t in tokens if t != "TER"]
    for idx, a in enumerate(new_atoms, start=1):
        a["serial"] = idx

    # CRYST1 box
    xs = [a["x"] for a in new_atoms]; ys = [a["y"] for a in new_atoms]; zs = [a["z"] for a in new_atoms]
    xlen = max(xs) - min(xs) + 2.0
    ylen = max(ys) - min(ys) + 2.0
    zlen = max(zs) - min(zs) + 2.0
    cryst1 = (f"CRYST1{xlen:9.3f}{ylen:9.3f}{zlen:9.3f}"
              f"  90.00  90.00  90.00 P 1           1\n")

    # CONECT for POPC residues
    residue_atoms = defaultdict(dict)  # (resname,chain,resnum,segid) -> {atomname: serial}
    for a in new_atoms:
        if a["resname"] == "POPC":
            key = (a["resname"], a["chain"], a["resnum"], a["segid"])
            residue_atoms[key][a["name"]] = a["serial"]

    all_conect = []
    n_popc = 0
    for key, sbn in residue_atoms.items():
        all_conect += make_conect(sbn, popc_atoms, popc_bonds)
        n_popc += 1

    # CONECT for the LEU46+cap junction bonds (not covered by OpenMM's
    # name-based standard-bonds builder, since NAM/HAM/CAM/HAM1-3 aren't
    # recognized names within its built-in "LEU" bonding table).
    cap_pairs = []
    for c_atom, nam, ham, cam, ham1, ham2, ham3 in cap_junctions:
        cap_pairs += [(c_atom["serial"], nam["serial"]),
                      (nam["serial"], ham["serial"]),
                      (nam["serial"], cam["serial"]),
                      (cam["serial"], ham1["serial"]),
                      (cam["serial"], ham2["serial"]),
                      (cam["serial"], ham3["serial"])]
    all_conect += conect_from_serial_pairs(cap_pairs)

    with open(out_path, "w") as f:
        f.write(cryst1)
        for t in tokens:
            f.write("TER\n" if t == "TER" else fmt_atom(t))
        for cl in all_conect:
            f.write(cl)
        f.write("END\n")

    print(f"         atoms       : {len(new_atoms)} atoms, {sum(1 for t in tokens if t == 'TER')} TER records")
    print(f"         renames     : " + "  ".join(f"{k}->{v}: {renamed_count[k]}" for k, v in RENAME.items()))
    print(f"         caps fixed  : {n_caps_fixed} chains (-CONH2 -> N-methylamide, fused into LEU46)")
    print(f"         box         : {xlen:.2f} x {ylen:.2f} x {zlen:.2f} A")
    print(f"         CONECT      : {n_popc} POPC -> {len(all_conect)} lines")
    print(f"         output      : {out_path}")
    print()

print("Done.")
