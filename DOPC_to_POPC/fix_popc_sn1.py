"""
Fix POPC sn-1 palmitoyl chain multipole frames.
Creates custom types 10004001-10004018 from POPE's palmitoyl parameters
with corrected frame references for the DOPC/POPE junction at C37-C38.
"""
import re

def read_file(path):
    with open(path) as f:
        return f.read()

pope_xml = read_file('/tab3/aman/AMEOBA_project/popc_from_published/pope.xml')
popc_xml = read_file('/tab3/aman/AMEOBA_project/popc_from_published/popc.xml')

# ── Type mapping ─────────────────────────────────────────────────────────────
pope_types = list(range(10002736, 10002754))   # 18 pope types for C38-C316+H's
custom_types = list(range(10004001, 10004019))
p2c = dict(zip(pope_types, custom_types))       # pope → custom
c2p = dict(zip(custom_types, pope_types))       # custom → pope

def translate_frame_val(pope_type, attr_name, val):
    """Map frame kx/kz from POPE type space to custom type space.
    Special: kz of C38 (10002736) → 1000341 (DOPC C37 type)."""
    ival = int(val)
    if pope_type == 10002736 and attr_name == 'kz' and ival == 10002734:
        return '1000341'   # kz of C38 → DOPC C37
    if ival in p2c:
        return str(p2c[ival])
    return val

# ── Extract POPE sections ─────────────────────────────────────────────────────
def get_section(xml, tag):
    m = re.search(rf'<{tag}[^>]*>(.*?)</{tag}>', xml, re.DOTALL)
    return m.group(0) if m else ''

pope_atomtypes  = get_section(pope_xml, 'AtomTypes')
pope_multipole  = get_section(pope_xml, 'AmoebaMultipoleForce')
pope_polarize   = get_section(pope_xml, 'AmoebaPolarizeForce')
pope_vdw        = get_section(pope_xml, 'AmoebaVdwForce')
pope_bond       = get_section(pope_xml, 'AmoebaBondForce')
pope_angle      = get_section(pope_xml, 'AmoebaAngleForce')
pope_torsion    = get_section(pope_xml, 'AmoebaTorsionForce')

# ── Build new AtomType entries ────────────────────────────────────────────────
new_atomtypes = []
for pt in pope_types:
    ct = p2c[pt]
    m = re.search(rf'<Type\s+name="{pt}"[^>]*/>', pope_atomtypes)
    if m:
        line = m.group().replace(f'name="{pt}"', f'name="{ct}"')
        new_atomtypes.append('   ' + line)

# ── Build new Multipole entries ───────────────────────────────────────────────
multipole_entries = {}
for m in re.finditer(r'<Multipole\s+[^>]*/>', pope_multipole):
    line = m.group()
    t = re.search(r'type="(\d+)"', line)
    if t and int(t.group(1)) in p2c:
        multipole_entries[int(t.group(1))] = line

new_multipoles = []
for pt in pope_types:
    ct = p2c[pt]
    line = multipole_entries[pt]
    # Replace type
    line = re.sub(r'type="\d+"', f'type="{ct}"', line)
    # Replace kx/kz frame references
    def replace_frame(m2):
        attr = m2.group(1)
        val  = m2.group(2)
        new_val = translate_frame_val(pt, attr, val)
        return f'{attr}="{new_val}"'
    line = re.sub(r'(k[xz])="(-?\d+)"', replace_frame, line)
    new_multipoles.append('   ' + line)

# ── Extract POPE VdW entries for the custom types ────────────────────────────
new_vdw = []
for m in re.finditer(r'<Vdw\s+[^>]*/>', pope_vdw):
    line = m.group()
    t = re.search(r'type="(\d+)"', line)
    if t and int(t.group(1)) in p2c:
        pt = int(t.group(1))
        new_line = line.replace(f'type="{pt}"', f'type="{p2c[pt]}"')
        new_vdw.append('   ' + new_line)

# ── Extract POPE Polarize entries ─────────────────────────────────────────────
new_polarize = []
for m in re.finditer(r'<Polarize\s+[^>]*/>', pope_polarize):
    line = m.group()
    t = re.search(r'type="(\d+)"', line)
    if t and int(t.group(1)) in p2c:
        pt = int(t.group(1))
        new_line = line.replace(f'type="{pt}"', f'type="{p2c[pt]}"')
        # Also replace any type references inside (polarize group)
        for p, c in p2c.items():
            new_line = new_line.replace(f'"{p}"', f'"{c}"')
        new_polarize.append('   ' + new_line)

# ── Extract POPE Bond/Angle/Torsion entries for classes 10002719,10002720,10002721,10002722 ──
pope_classes = {'10002719', '10002720', '10002721', '10002722'}

def extract_force_entries(section, tag, n_classes):
    """Extract entries where ANY class attr is in pope_classes."""
    entries = []
    pattern = rf'<{tag}\s+[^>]*/>'
    for m in re.finditer(pattern, section):
        line = m.group()
        classes_in_line = set(re.findall(r'class\d*="(\d+)"', line))
        if classes_in_line & pope_classes:
            entries.append('   ' + line)
    return entries

new_bonds    = extract_force_entries(pope_bond,    'Bond',    2)
new_angles   = extract_force_entries(pope_angle,   'Angle',   3)
new_torsions = extract_force_entries(pope_torsion, 'Torsion', 4)

print(f"New entries: {len(new_atomtypes)} types, {len(new_multipoles)} multipoles, "
      f"{len(new_vdw)} vdw, {len(new_polarize)} polarize, "
      f"{len(new_bonds)} bonds, {len(new_angles)} angles, {len(new_torsions)} torsions")

# ── Cross-junction parameters: DOPC C37 (class 1000323) ↔ C38 (class 10002719) ──
# Get POPE's 10002719-10002719 C-C bond params for cross-junction
cross_bond = None
for m in re.finditer(r'<Bond\s+[^>]*/>', pope_bond):
    line = m.group()
    c1 = re.search(r'class1="(\d+)"', line)
    c2 = re.search(r'class2="(\d+)"', line)
    if c1 and c2 and c1.group(1)=='10002719' and c2.group(1)=='10002719':
        cross_bond = line.replace('class1="10002719"', 'class1="1000323"')
        cross_bond = cross_bond.replace('class2="10002719"', 'class2="10002719"')
        break

# Cross-junction angles
cross_angles = []
for m in re.finditer(r'<Angle\s+[^>]*/>', pope_angle):
    line = m.group()
    cs = re.findall(r'class\d+="(\d+)"', line)
    if all(c == '10002719' for c in cs):
        # sp3 C-C-C angle, make cross version: 1000323-10002719-10002719
        a1 = line.replace('class1="10002719"', 'class1="1000323"', 1)
        cross_angles.append('   ' + a1)
        # and 10002719-1000323-10002719 ... no, that doesn't exist in this system
    # Also need H-C-C angles: 10002720-10002719-10002719
    if cs == ['10002720','10002719','10002719'] or cs == ['10002719','10002719','10002720']:
        # 1000324-1000323-10002719 (H7-C37-C38) and 1000323-10002719-10002720 (C37-C38-H8)
        a1 = line.replace('class1="10002720"','class1="1000324"',1).replace('class2="10002719"','class2="1000323"',1)
        cross_angles.append('   ' + a1)
        a2 = line.replace('class1="10002719"','class1="1000323"',1)
        cross_angles.append('   ' + a2)

# Cross-junction torsions: use POPE's generic 10002719×4 torsion for all cross cases
pope_sp3_torsion = None
for m in re.finditer(r'<Torsion\s+[^>]*/>', pope_torsion):
    line = m.group()
    cs = re.findall(r'class\d+="(\d+)"', line)
    if all(c == '10002719' for c in cs):
        pope_sp3_torsion = line
        break
if pope_sp3_torsion is None:
    # Fallback: look for any 10002719-10002719-10002719-10002719
    print("WARNING: Could not find 10002719×4 torsion, using placeholder")
    pope_sp3_torsion = '<Torsion class1="10002719" class2="10002719" class3="10002719" class4="10002719" k1="0" k2="0" k3="0.3347"/>'

def make_cross_torsion(c1, c2, c3, c4):
    return pope_sp3_torsion.replace('class1="10002719"', f'class1="{c1}"',1)\
                            .replace('class2="10002719"', f'class2="{c2}"',1)\
                            .replace('class3="10002719"', f'class3="{c3}"',1)\
                            .replace('class4="10002719"', f'class4="{c4}"',1)

cross_torsions = [
    '   ' + make_cross_torsion('1000323','1000323','10002719','10002719'),  # C36-C37-C38-C39
    '   ' + make_cross_torsion('1000323','10002719','10002719','10002719'), # C37-C38-C39-C310
    '   ' + make_cross_torsion('1000324','1000323','10002719','10002719'),  # H7-C37-C38-C39
    '   ' + make_cross_torsion('1000323','1000323','10002719','10002720'),  # C36-C37-C38-H8
    '   ' + make_cross_torsion('1000324','1000323','10002719','10002720'),  # H7-C37-C38-H8
    '   ' + make_cross_torsion('1000323','10002719','10002719','10002720'), # C37-C38-C39-H9
]

print(f"Cross-junction: 1 bond, {len(cross_angles)} angles, {len(cross_torsions)} torsions")

# ── Atom name → new custom type mapping ──────────────────────────────────────
# POPE atom → POPE type → custom type
pope_res = re.search(r'<Residue name="POPE">(.*?)</Residue>', pope_xml, re.DOTALL).group(1)
atom_to_pope = {}
for m in re.finditer(r'<Atom name="([^"]+)" type="(\d+)"', pope_res):
    atom_to_pope[m.group(1)] = int(m.group(2))

# Atoms in sn-1 palmitoyl chain (C38-C316 and H's)
sn1_atoms_pope_types = {}
for atom_name, pope_type in atom_to_pope.items():
    if pope_type in p2c:
        sn1_atoms_pope_types[atom_name] = p2c[pope_type]

print("\nAtom → custom type in sn-1 palmitoyl chain:")
for an, ct in sorted(sn1_atoms_pope_types.items()):
    print(f"  {an}: {ct}")

# ── Patch POPC residue in popc.xml ───────────────────────────────────────────
# Update atom type assignments for sn-1 C38-C316 and H's
def update_popc_residue(xml, atom_to_new_type):
    """Update atom type assignments in POPC residue."""
    popc_res_m = re.search(r'(<Residue name="POPC">)(.*?)(</Residue>)', xml, re.DOTALL)
    if not popc_res_m:
        return xml, False
    pre, res_body, post = popc_res_m.group(1), popc_res_m.group(2), popc_res_m.group(3)
    changed = 0
    for atom_name, new_type in atom_to_new_type.items():
        # Match Atom name="C38" type="..."
        old_m = re.search(rf'<Atom name="{re.escape(atom_name)}" type="(\d+)"', res_body)
        if old_m:
            old_type = old_m.group(1)
            res_body = res_body.replace(
                f'<Atom name="{atom_name}" type="{old_type}"',
                f'<Atom name="{atom_name}" type="{new_type}"'
            )
            changed += 1
        else:
            print(f"  WARNING: atom {atom_name} not found in POPC residue")
    print(f"  Updated {changed} atom types in POPC residue")
    return xml.replace(popc_res_m.group(0), pre + res_body + post), changed > 0

# ── Patch force field sections ────────────────────────────────────────────────
def inject_before_close(xml, section_tag, new_lines):
    """Inject lines before closing tag of a section."""
    close_tag = f'</{section_tag}>'
    if close_tag not in xml:
        return xml
    inject = '\n'.join(new_lines) + '\n'
    return xml.replace(close_tag, inject + close_tag, 1)

print("\nPatching popc.xml...")

# 1. Add AtomTypes
popc_xml = inject_before_close(popc_xml, 'AtomTypes', new_atomtypes)

# 2. Add Multipoles
popc_xml = inject_before_close(popc_xml, 'AmoebaMultipoleForce', new_multipoles)

# 3. Add VdW
popc_xml = inject_before_close(popc_xml, 'AmoebaVdwForce', new_vdw)

# 4. Add Polarize
popc_xml = inject_before_close(popc_xml, 'AmoebaPolarizeForce', new_polarize)

# 5. Add Bonds (including cross-junction)
all_bonds = new_bonds + ([('   ' + cross_bond)] if cross_bond else [])
popc_xml = inject_before_close(popc_xml, 'AmoebaBondForce', all_bonds)

# 6. Add Angles (including cross-junction)
all_angles = new_angles + cross_angles
popc_xml = inject_before_close(popc_xml, 'AmoebaAngleForce', all_angles)

# 7. Add Torsions (including cross-junction)
all_torsions = new_torsions + cross_torsions
popc_xml = inject_before_close(popc_xml, 'AmoebaTorsionForce', all_torsions)

# 8. Update POPC residue atom types
popc_xml, ok = update_popc_residue(popc_xml, sn1_atoms_pope_types)

# ── Write output ──────────────────────────────────────────────────────────────
out_path = '/tab3/aman/AMEOBA_project/popc_from_published/popc.xml'
with open(out_path, 'w') as f:
    f.write(popc_xml)
print(f"\nWrote {out_path}")

