"""
Fix POPC sn-1 palmitoyl chain multipole frames (v2).
Creates custom types 10004001-10004018 from POPE's palmitoyl parameters
with corrected frame references at the DOPC/POPE junction (C37-C38).
"""
import re

def read_file(path):
    with open(path) as f:
        return f.read()

# Read ORIGINAL pope.xml and popc.xml (popc may have been partially patched)
pope_xml = read_file('/tab3/aman/AMEOBA_project/popc_from_published/pope.xml')
# Read from backup or revert to clean version
popc_xml = read_file('/tab3/aman/AMEOBA_project/popc_from_published/popc_backup.xml')

# ── Type mapping ─────────────────────────────────────────────────────────────
pope_types  = list(range(10002736, 10002754))
custom_types= list(range(10004001, 10004019))
p2c = dict(zip(pope_types, custom_types))
c2p = dict(zip(custom_types, pope_types))

def translate_frame_val(pope_type, attr_name, val):
    ival = int(val)
    # Special: C38 kz → DOPC C37 (type 1000341)
    if pope_type == 10002736 and attr_name == 'kz' and ival == 10002734:
        return '1000341'
    if ival in p2c:
        return str(p2c[ival])
    return val

# ── Helper: extract section (handles multi-attribute opening tags) ─────────────
def get_section(xml, tag):
    m = re.search(rf'<{tag}[\s>][^>]*>(.*?)</{tag}>', xml, re.DOTALL)
    if not m:
        m = re.search(rf'<{tag}>(.*?)</{tag}>', xml, re.DOTALL)
    return m.group(0) if m else ''

pope_atomtypes = get_section(pope_xml, 'AtomTypes')
pope_multipole = get_section(pope_xml, 'AmoebaMultipoleForce')
pope_vdw       = get_section(pope_xml, 'AmoebaVdwForce')
pope_bond      = get_section(pope_xml, 'AmoebaBondForce')
pope_angle     = get_section(pope_xml, 'AmoebaAngleForce')
pope_torsion   = get_section(pope_xml, 'AmoebaTorsionForce')

# ── 1. New AtomType entries ───────────────────────────────────────────────────
# Type entries have name= as last attribute: <Type class="..." element="..." mass="..." name="10002736"/>
new_atomtypes = []
for pt in pope_types:
    ct = p2c[pt]
    # Search for name="${pt}" anywhere in the Type entry
    m = re.search(rf'<Type\s+[^>]*name="{pt}"[^>]*/>', pope_atomtypes)
    if m:
        line = m.group().replace(f'name="{pt}"', f'name="{ct}"')
        new_atomtypes.append('      ' + line)
    else:
        print(f"WARNING: type {pt} not found in AtomTypes")

print(f"AtomTypes: {len(new_atomtypes)}/18 found")

# ── 2. New Multipole entries ─────────────────────────────────────────────────
new_multipoles = []
for pt in pope_types:
    ct = p2c[pt]
    m = re.search(rf'<Multipole\s+[^>]*type="{pt}"[^>]*/>', pope_multipole)
    if m:
        line = m.group()
        line = re.sub(r'type="\d+"', f'type="{ct}"', line)
        def replace_frame(m2):
            attr, val = m2.group(1), m2.group(2)
            return f'{attr}="{translate_frame_val(pt, attr, val)}"'
        line = re.sub(r'(k[xz])="(-?\d+)"', replace_frame, line)
        new_multipoles.append('      ' + line)
    else:
        print(f"WARNING: type {pt} not found in Multipoles")

print(f"Multipoles: {len(new_multipoles)}/18")

# ── 3. VdW entries for POPE classes 10002719-10002722 (use class= attr) ──────
pope_classes_needed = {'10002719','10002720','10002721','10002722'}
new_vdw = []
seen_vdw_classes = set()
for m in re.finditer(r'<Vdw\s+[^>]*/>', pope_vdw):
    line = m.group()
    cl = re.search(r'class="(\d+)"', line)
    if cl and cl.group(1) in pope_classes_needed and cl.group(1) not in seen_vdw_classes:
        new_vdw.append('      ' + line)
        seen_vdw_classes.add(cl.group(1))

print(f"VdW: {len(new_vdw)} entries")

# ── 4. Polarize entries (no pgrp — each atom in its own group, approximation) ─
# Carbon sp3 CH2: polarizability=0.001334, thole=0.390
# Hydrogen sp3 CH2-H: polarizability=0.000496, thole=0.390
# terminal CH3 C: polarizability=0.001334, H: 0.000496
carbon_types = {10004001,10004003,10004005,10004007,10004009,10004011,10004013,10004015,10004017}
new_polarize = []
for ct in custom_types:
    polz = '0.0013340000000000001' if ct in carbon_types else '0.000496'
    new_polarize.append(f'      <Polarize polarizability="{polz}" thole="0.390" type="{ct}"/>')

print(f"Polarize: {len(new_polarize)} entries")

# ── 5. Bond entries for POPE classes ─────────────────────────────────────────
def extract_by_class(section, tag, classes_needed):
    entries = []
    seen = set()
    for m in re.finditer(rf'<{tag}\s+[^>]*/>', section):
        line = m.group()
        cls = set(re.findall(r'class\d*="(\d+)"', line))
        if cls & classes_needed and line not in seen:
            entries.append('      ' + line)
            seen.add(line)
    return entries

new_bonds   = extract_by_class(pope_bond,   'Bond',  pope_classes_needed)
new_angles  = extract_by_class(pope_angle,  'Angle', pope_classes_needed)
# Torsion elements inside AmoebaTorsionForce are named <Angle>!
new_torsions = extract_by_class(pope_torsion, 'Angle', pope_classes_needed)

print(f"Bond: {len(new_bonds)}, Angle: {len(new_angles)}, Torsion(Angle): {len(new_torsions)}")

# ── 6. Cross-junction parameters: DOPC C37 (class 1000323) ↔ C38 (class 10002719) ──
# Bond: 1000323-10002719
cross_bond = None
for m in re.finditer(r'<Bond\s+[^>]*/>', pope_bond):
    line = m.group()
    c1 = re.search(r'class1="(\d+)"', line)
    c2 = re.search(r'class2="(\d+)"', line)
    if c1 and c2 and c1.group(1)=='10002719' and c2.group(1)=='10002719':
        cross_bond = '      ' + line.replace('class1="10002719"','class1="1000323"',1)
        break

# Cross angles (use POPE 10002719-10002719-10002719 angle for all sp3 cross-angles)
pope_sp3_angle = None
for m in re.finditer(r'<Angle\s+[^>]*/>', pope_angle):
    line = m.group()
    cs = re.findall(r'class\d+="(\d+)"', line)
    if all(c=='10002719' for c in cs):
        pope_sp3_angle = line
        break
pope_cch_angle = None
for m in re.finditer(r'<Angle\s+[^>]*/>', pope_angle):
    line = m.group()
    cs = re.findall(r'class\d+="(\d+)"', line)
    if cs == ['10002720','10002719','10002719']:
        pope_cch_angle = line
        break

cross_angles = []
if pope_sp3_angle:
    # C36-C37-C38: class3 10002719 → class3 stays, class1=1000323, class2=1000323
    a = pope_sp3_angle.replace('class1="10002719"','class1="1000323"',1).replace('class2="10002719"','class2="1000323"',1)
    cross_angles.append('      ' + a)
    # C37-C38-C39: class1=1000323, class2=10002719, class3=10002719
    a = pope_sp3_angle.replace('class1="10002719"','class1="1000323"',1)
    cross_angles.append('      ' + a)
if pope_cch_angle:
    # H7-C37-C38: 1000324-1000323-10002719
    a = pope_cch_angle.replace('class2="10002719"','class2="1000323"',1)
    cross_angles.append('      ' + a)
    # C37-C38-H8: 1000323-10002719-10002720
    a = pope_cch_angle.replace('class1="10002720"','class1="1000323"',1).replace('class2="10002719" ','class2="10002719" ',1)
    # This swaps class1/class3 but the angle is symmetric
    a2 = pope_cch_angle.replace('class1="10002720"','class1="1000323"',1)
    cross_angles.append('      ' + a2)

# Cross torsions — use POPE sp3 4×10002719 torsion
pope_sp3_torsion = None
for m in re.finditer(r'<Angle\s+[^>]*/>', pope_torsion):
    line = m.group()
    cs = re.findall(r'class\d+="(\d+)"', line)
    if all(c=='10002719' for c in cs):
        pope_sp3_torsion = line
        break

cross_torsions = []
if pope_sp3_torsion:
    def xtor(c1,c2,c3,c4):
        return '      ' + pope_sp3_torsion\
            .replace('class1="10002719"',f'class1="{c1}"',1)\
            .replace('class2="10002719"',f'class2="{c2}"',1)\
            .replace('class3="10002719"',f'class3="{c3}"',1)\
            .replace('class4="10002719"',f'class4="{c4}"',1)
    cross_torsions = [
        xtor('1000323','1000323','10002719','10002719'),  # C36-C37-C38-C39
        xtor('1000323','10002719','10002719','10002719'), # C37-C38-C39-C310
        xtor('1000324','1000323','10002719','10002719'),  # H7-C37-C38-C39
        xtor('1000323','1000323','10002719','10002720'),  # C36-C37-C38-H8
        xtor('1000324','1000323','10002719','10002720'),  # H7-C37-C38-H8
        xtor('1000323','10002719','10002719','10002720'), # C37-C38-C39-H9
    ]

print(f"Cross: 1 bond, {len(cross_angles)} angles, {len(cross_torsions)} torsions")

# ── 7. Get atom→custom_type mapping from POPE residue ────────────────────────
pope_res = re.search(r'<Residue name="POPE">(.*?)</Residue>', pope_xml, re.DOTALL).group(1)
atom_to_new_type = {}
for m in re.finditer(r'<Atom name="([^"]+)" type="(\d+)"', pope_res):
    aname, pt = m.group(1), int(m.group(2))
    if pt in p2c:
        atom_to_new_type[aname] = p2c[pt]

print(f"\nAtoms to retype: {len(atom_to_new_type)}")
for an, ct in sorted(atom_to_new_type.items()):
    print(f"  {an}: pope {c2p[ct]} → custom {ct}")

# ── 8. Inject into popc.xml ───────────────────────────────────────────────────
def inject_before_close(xml, tag, lines):
    close_tag = f'</{tag}>'
    if close_tag not in xml:
        # Try with leading whitespace
        for line in xml.split('\n'):
            if close_tag in line:
                close_tag = line.strip()  # won't work well
        print(f"WARNING: close tag </{tag}> not found directly")
        return xml
    insert = '\n'.join(lines) + '\n'
    return xml.replace(close_tag, insert + close_tag, 1)

def update_popc_residue(xml, atom_to_new_type):
    popc_res_m = re.search(r'(<Residue name="POPC">)(.*?)(</Residue>)', xml, re.DOTALL)
    if not popc_res_m:
        print("WARNING: POPC residue not found")
        return xml, 0
    pre, res_body, post = popc_res_m.group(1), popc_res_m.group(2), popc_res_m.group(3)
    changed = 0
    for aname, new_type in atom_to_new_type.items():
        old_m = re.search(rf'<Atom name="{re.escape(aname)}" type="(\d+)"', res_body)
        if old_m:
            old_type = old_m.group(1)
            res_body = res_body.replace(
                f'<Atom name="{aname}" type="{old_type}"',
                f'<Atom name="{aname}" type="{new_type}"'
            )
            changed += 1
        else:
            print(f"  WARN: atom {aname} not in POPC residue")
    return xml.replace(popc_res_m.group(0), pre + res_body + post), changed

print("\nPatching popc.xml...")
xml = popc_xml

xml = inject_before_close(xml, 'AtomTypes', new_atomtypes)
xml = inject_before_close(xml, 'AmoebaMultipoleForce', new_multipoles)
xml = inject_before_close(xml, 'AmoebaVdwForce', new_vdw)
xml = inject_before_close(xml, 'AmoebaPolarizeForce', new_polarize)
xml = inject_before_close(xml, 'AmoebaBondForce', new_bonds + ([cross_bond] if cross_bond else []))
xml = inject_before_close(xml, 'AmoebaAngleForce', new_angles + cross_angles)
xml = inject_before_close(xml, 'AmoebaTorsionForce', new_torsions + cross_torsions)

xml, n_changed = update_popc_residue(xml, atom_to_new_type)
print(f"Updated {n_changed} atom types in POPC residue")

# Write output
out_path = '/tab3/aman/AMEOBA_project/popc_from_published/popc.xml'
with open(out_path, 'w') as f:
    f.write(xml)
print(f"\nWrote {out_path}")

