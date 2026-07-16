"""
Patch popc.xml to add missing pieces:
1. AtomType entries for custom types 10004001-10004018
2. VdW entries for POPE classes 10002719, 10002720, 10002721, 10002722
3. Polarize entries for custom types 10004001-10004018
4. Torsion (Angle-tagged) entries for POPE sp3 chain + cross-junction
"""
import re

def read_file(path):
    with open(path) as f:
        return f.read()

pope_xml = read_file('/tab3/aman/AMEOBA_project/popc_from_published/pope.xml')
popc_xml = read_file('/tab3/aman/AMEOBA_project/popc_from_published/popc.xml')

pope_types   = list(range(10002736, 10002754))
custom_types = list(range(10004001, 10004019))
p2c = dict(zip(pope_types, custom_types))

carbon_custom = {10004001,10004003,10004005,10004007,10004009,10004011,10004013,10004015,10004017}

def get_section_inner(xml, tag):
    m = re.search(rf'<{tag}[\s>][^>]*>(.*?)</{tag}>', xml, re.DOTALL)
    if not m:
        m = re.search(rf'<{tag}>(.*?)</{tag}>', xml, re.DOTALL)
    return m.group(1) if m else ''

# ── 1. AtomType entries ────────────────────────────────────────────────────────
pope_at_inner = get_section_inner(pope_xml, 'AtomTypes')
new_atomtypes = []
for pt in pope_types:
    ct = p2c[pt]
    # name= is last attribute: <Type class="..." element="..." mass="..." name="10002736"/>
    m = re.search(rf'<Type\s+[^>]*\bname="{pt}"[^>]*/>', pope_at_inner)
    if m:
        line = m.group().replace(f'name="{pt}"', f'name="{ct}"')
        new_atomtypes.append('      ' + line)
    else:
        print(f"WARNING: AtomType {pt} not found")
print(f"AtomTypes: {len(new_atomtypes)}/18")

# ── 2. VdW entries for POPE classes 10002719,10002720,10002721,10002722 ───────
pope_vdw_inner = get_section_inner(pope_xml, 'AmoebaVdwForce')
pope_classes_needed = {'10002719','10002720','10002721','10002722'}
new_vdw = []
for m in re.finditer(r'<Vdw\s+[^>]*/>', pope_vdw_inner):
    line = m.group()
    cl = re.search(r'class="(\d+)"', line)
    if cl and cl.group(1) in pope_classes_needed:
        new_vdw.append('      ' + line)
print(f"VdW: {len(new_vdw)} entries")

# ── 3. Polarize entries ───────────────────────────────────────────────────────
new_polarize = []
for ct in custom_types:
    polz = '0.0013340000000000001' if ct in carbon_custom else '0.000496'
    new_polarize.append(f'      <Polarize polarizability="{polz}" thole="0.390" type="{ct}"/>')
print(f"Polarize: {len(new_polarize)}")

# ── 4. Torsion entries (tag is <Angle> inside AmoebaTorsionForce) ─────────────
pope_tor_inner = get_section_inner(pope_xml, 'AmoebaTorsionForce')
new_torsions = []
seen = set()
for m in re.finditer(r'<Angle\s+[^>]*/>', pope_tor_inner):
    line = m.group()
    cls = set(re.findall(r'class\d+="(\d+)"', line))
    if cls & pope_classes_needed and line not in seen:
        new_torsions.append('      ' + line)
        seen.add(line)
print(f"Torsions: {len(new_torsions)}")

# Cross-junction torsions (C37 DOPC class 1000323 ↔ C38 POPE class 10002719)
pope_sp3_tor = None
for m in re.finditer(r'<Angle\s+[^>]*/>', pope_tor_inner):
    line = m.group()
    cs = re.findall(r'class\d+="(\d+)"', line)
    if all(c=='10002719' for c in cs):
        pope_sp3_tor = line
        break

cross_tors = []
if pope_sp3_tor:
    def xtor(c1,c2,c3,c4):
        return '      ' + pope_sp3_tor\
            .replace('class1="10002719"',f'class1="{c1}"',1)\
            .replace('class2="10002719"',f'class2="{c2}"',1)\
            .replace('class3="10002719"',f'class3="{c3}"',1)\
            .replace('class4="10002719"',f'class4="{c4}"',1)
    cross_tors = [
        xtor('1000323','1000323','10002719','10002719'),
        xtor('1000323','10002719','10002719','10002719'),
        xtor('1000324','1000323','10002719','10002719'),
        xtor('1000323','1000323','10002719','10002720'),
        xtor('1000324','1000323','10002719','10002720'),
        xtor('1000323','10002719','10002719','10002720'),
    ]
    print(f"Cross torsions: {len(cross_tors)}")
else:
    print("WARNING: no 10002719×4 torsion found")

# ── Inject into popc.xml ───────────────────────────────────────────────────────
def inject_before_close(xml, tag, lines):
    close_tag = f'</{tag}>'
    insert = '\n'.join(lines) + '\n'
    return xml.replace(close_tag, insert + close_tag, 1)

xml = popc_xml
xml = inject_before_close(xml, 'AtomTypes', new_atomtypes)
xml = inject_before_close(xml, 'AmoebaVdwForce', new_vdw)
xml = inject_before_close(xml, 'AmoebaPolarizeForce', new_polarize)
xml = inject_before_close(xml, 'AmoebaTorsionForce', new_torsions + cross_tors)

out_path = '/tab3/aman/AMEOBA_project/popc_from_published/popc.xml'
with open(out_path, 'w') as f:
    f.write(xml)
print(f"\nWrote {out_path} ({len(xml.splitlines())} lines)")

