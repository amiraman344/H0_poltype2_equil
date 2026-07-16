"""
Compare poltype2-derived POPC force field (final.key) vs H0 DOPC-based popc.xml.

Focuses on key AMOEBA parameters: multipoles, polarizabilities, vdW.
"""
import re, math, xml.etree.ElementTree as ET

KEY_FILE = "/tab3/aman/AMEOBA_project/poltype2_popc/final.key"
XML_FILE = "/tab1/aaamir2/Christian_Openmm_project/AMOEBA/openmm_simulation/H0/popc.xml"

# ─── parse poltype2 final.key ────────────────────────────────────────────────

def parse_key(path):
    atoms      = {}  # type -> {sym, element, mass, valence}
    polarize   = {}  # type -> {alpha, thole}
    vdw        = {}  # type -> {rmin, eps}
    multipoles = {}  # type -> {charge, dipole[3], quad[6]}

    with open(path) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("#"):
            i += 1; continue

        tok = line.split()
        kw = tok[0].lower()

        if kw == "atom":
            # atom  401  401  N  "desc"  7  14.007  4
            m = re.match(r'atom\s+(\d+)\s+\d+\s+(\S+)\s+"[^"]*"\s+(\d+)\s+([\d.]+)\s+(\d+)', line, re.I)
            if m:
                t = int(m.group(1))
                atoms[t] = {"sym": m.group(2), "element": m.group(3),
                             "mass": float(m.group(4)), "val": int(m.group(5))}

        elif kw == "polarize":
            t = int(tok[1])
            polarize[t] = {"alpha": float(tok[2]), "thole": float(tok[3])}

        elif kw == "vdw":
            t = int(tok[1])
            vdw[t] = {"rmin": float(tok[2]), "eps": float(tok[3])}

        elif kw == "multipole":
            t = int(tok[1])
            charge = float(tok[-1])
            if i+4 < len(lines):
                d = list(map(float, lines[i+1].split()))
                q1 = list(map(float, lines[i+2].split()))
                q2 = list(map(float, lines[i+3].split()))
                q3 = list(map(float, lines[i+4].split()))
                quad = q1 + q2 + q3
            else:
                d = [0,0,0]; quad = [0]*6
            multipoles[t] = {"charge": charge, "dipole": d, "quad": quad}
            i += 4

        i += 1

    return atoms, polarize, vdw, multipoles


# ─── parse H0 popc.xml ───────────────────────────────────────────────────────

def parse_xml(path):
    tree = ET.parse(path)
    root = tree.getroot()

    atom_types = {}  # name -> {class, element, mass}
    polarize_xml = {}  # class -> {alpha_A3}
    vdw_xml     = {}  # class -> {sigma_A, eps_kcal}
    multipoles_xml = {}  # type-name -> {charge, dipole, quad}

    for at in root.iter("Type"):
        name = at.get("name"); cls = at.get("class")
        elem = at.get("element","?"); mass = float(at.get("mass", 0))
        atom_types[name] = {"class": cls, "element": elem, "mass": mass}

    KJ_TO_KCAL = 1.0 / 4.184
    NM_TO_A = 10.0

    for vf in root.iter("AmoebaVdwForce"):
        for v in vf.iter("Vdw"):
            cls = v.get("class")
            sigma_nm = float(v.get("sigma"))   # nm, diameter=Rmin in AMOEBA
            eps_kj   = float(v.get("epsilon")) # kJ/mol
            vdw_xml[cls] = {
                "rmin_A": sigma_nm * NM_TO_A,          # Rmin in Angstrom
                "eps_kcal": eps_kj * KJ_TO_KCAL
            }

    for pf in root.iter("AmoebaMultipoleForce"):
        for mp in pf.iter("Multipole"):
            tname = mp.get("type")
            c  = float(mp.get("c0"))
            dx = float(mp.get("d1")); dy = float(mp.get("d2")); dz = float(mp.get("d3"))
            q = [float(mp.get(f"q{j}",0)) for j in ["11","21","22","31","32","33"]]
            multipoles_xml[tname] = {"charge": c, "dipole": [dx,dy,dz], "quad": q}

    for pf in root.iter("AmoebaGeneralizedKirkwoodForce"):
        pass  # not needed here

    # polarize — lives inside AmoebaMultipoleForce in this XML format
    for p in root.iter("Polarize"):
        tname = p.get("type")
        alpha_nm3 = float(p.get("polarizability"))  # nm^3 → Å^3 × 1000
        thole = float(p.get("thole", 0.39))
        cls = atom_types.get(tname, {}).get("class", tname)
        polarize_xml[tname] = {"alpha": alpha_nm3 * 1000, "thole": thole, "cls": cls}

    return atom_types, polarize_xml, vdw_xml, multipoles_xml


# ─── aggregate by element ────────────────────────────────────────────────────

def stats(vals):
    if not vals: return 0, 0, 0
    mn, mx = min(vals), max(vals)
    avg = sum(vals)/len(vals)
    return avg, mn, mx


def by_element(atoms, param_dict, key):
    from collections import defaultdict
    el_vals = defaultdict(list)
    for t, info in atoms.items():
        el = info["element"]
        v = param_dict.get(t, {}).get(key)
        if v is not None:
            el_vals[el].append(v)
    return el_vals


# ─── main comparison ─────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("POPC AMOEBA Force Field Comparison")
    print("  poltype2 (QM-derived)  vs  H0 DOPC-based (dopc.prm)")
    print("=" * 70)

    atoms_p, pol_p, vdw_p, mp_p = parse_key(KEY_FILE)
    atom_types_x, pol_x, vdw_x, mp_x = parse_xml(XML_FILE)

    # ── atom count ──
    print(f"\n[Atom Types]")
    print(f"  poltype2 : {len(atoms_p)} unique types ({len(atoms_p)} atoms in POPC)")
    print(f"  H0 xml   : {len(atom_types_x)} unique types (shared DOPC classes)")

    # Atomic number → symbol mapping
    ANUM = {"1":"H","6":"C","7":"N","8":"O","15":"P","16":"S","9":"F"}
    for t in atoms_p:
        atoms_p[t]["element"] = ANUM.get(atoms_p[t]["element"], atoms_p[t]["element"])

    el_count_p = {}
    for t, info in atoms_p.items():
        el_count_p[info["element"]] = el_count_p.get(info["element"], 0) + 1
    el_count_x = {}
    for n, info in atom_types_x.items():
        el = info["element"]
        el_count_x[el] = el_count_x.get(el, 0) + 1
    print(f"  poltype2 elements: {el_count_p}")
    print(f"  H0 xml elements  : {el_count_x}")

    # ── polarizability by element ──
    print(f"\n[Polarizability α (Å³)] — poltype2 vs H0")
    print(f"  {'Element':6}  {'P2 avg':>8} {'P2 min-max':>16}   {'H0 avg':>8} {'H0 min-max':>16}")
    print(f"  {'-'*64}")

    el_pol_p = {}
    for t, info in pol_p.items():
        el = atoms_p.get(t, {}).get("element", "?")
        el_pol_p.setdefault(el, []).append(info["alpha"])

    el_pol_x = {}
    for tname, info in pol_x.items():
        el = atom_types_x.get(tname, {}).get("element", "?")
        el_pol_x.setdefault(el, []).append(info["alpha"])

    for el in sorted(set(list(el_pol_p.keys()) + list(el_pol_x.keys()))):
        avg_p, mn_p, mx_p = stats(el_pol_p.get(el, []))
        avg_x, mn_x, mx_x = stats(el_pol_x.get(el, []))
        print(f"  {el:6}  {avg_p:8.4f} {mn_p:7.4f}-{mx_p:.4f}   {avg_x:8.4f} {mn_x:7.4f}-{mx_x:.4f}")

    # ── vdW by element ──
    print(f"\n[vdW] Rmin (Å) and ε (kcal/mol) — poltype2 vs H0")
    print(f"  {'Element':6}  {'P2 Rmin':>8} {'P2 eps':>8}   {'H0 Rmin':>8} {'H0 eps':>8}")
    print(f"  {'-'*56}")

    el_rmin_p = {}; el_eps_p = {}
    for t, v in vdw_p.items():
        el = atoms_p.get(t, {}).get("element", "?")
        el_rmin_p.setdefault(el, []).append(v["rmin"])
        el_eps_p.setdefault(el, []).append(v["eps"])

    el_rmin_x = {}; el_eps_x = {}
    for cls, v in vdw_x.items():
        # map class to element via atom_types
        el = next((info["element"] for info in atom_types_x.values()
                   if info["class"] == cls), "?")
        el_rmin_x.setdefault(el, []).append(v["rmin_A"])
        el_eps_x.setdefault(el, []).append(v["eps_kcal"])

    for el in sorted(set(list(el_rmin_p.keys()) + list(el_rmin_x.keys()))):
        avg_rp, _, _ = stats(el_rmin_p.get(el, []))
        avg_ep, _, _ = stats(el_eps_p.get(el, []))
        avg_rx, _, _ = stats(el_rmin_x.get(el, []))
        avg_ex, _, _ = stats(el_eps_x.get(el, []))
        print(f"  {el:6}  {avg_rp:8.4f} {avg_ep:8.4f}   {avg_rx:8.4f} {avg_ex:8.4f}")

    # ── Phosphorus specifically ──
    print(f"\n[Phosphorus detailed comparison]")
    P_type_p = next(t for t, a in atoms_p.items() if a["element"] == "P")
    P_name_x = next(n for n, a in atom_types_x.items() if a["element"] == "P")

    vp = vdw_p.get(P_type_p, {})
    vx = vdw_x.get(atom_types_x[P_name_x]["class"], {})
    print(f"  vdW Rmin  : poltype2 = {vp.get('rmin','?'):.4f} Å   H0 = {vx.get('rmin_A','?'):.4f} Å")
    print(f"  vdW eps   : poltype2 = {vp.get('eps','?'):.4f} kcal/mol   H0 = {vx.get('eps_kcal','?'):.4f} kcal/mol")
    pp = pol_p.get(P_type_p, {})
    px = pol_x.get(P_name_x, pol_x.get(str(P_name_x), {}))
    alpha_p = pp.get('alpha', float('nan'))
    alpha_x = px.get('alpha', float('nan'))
    print(f"  Polar α   : poltype2 = {alpha_p:.4f} Å³   H0 = {alpha_x:.4f} Å³")
    mp_val_p = mp_p.get(P_type_p, {})
    mp_val_x = mp_x.get(P_name_x, {})
    chg_p = mp_val_p.get('charge', float('nan'))
    chg_x = mp_val_x.get('charge', float('nan'))
    print(f"  Monopole  : poltype2 = {chg_p:.5f} e   H0 = {chg_x:.5f} e")
    dip_p = mp_val_p.get("dipole", [0,0,0])
    dip_x = mp_val_x.get("dipole", [0,0,0])
    # H0 dipole is in e·nm, convert to e·Å: * 10
    dip_x_A = [d * 10 for d in dip_x]
    print(f"  Dipole(Å) : poltype2 = {[f'{d:.5f}' for d in dip_p]}   H0 = {[f'{d:.5f}' for d in dip_x_A]}")

    # ── charge balance ──
    print(f"\n[Total charge (monopole sum)]")
    total_p = sum(v["charge"] for v in mp_p.values())
    total_x = sum(v["charge"] for v in mp_x.values())
    print(f"  poltype2  : {total_p:.4f} e  (expected 0 for neutral POPC)")
    print(f"  H0 xml    : {total_x:.4f} e")

    # ── molecular dipole magnitude ──
    print(f"\n[Molecular dipole from permanent monopoles+dipoles (gas phase estimate)]")
    # Sum dipole contributions (local frame — approximate, not accounting for positions)
    # Just show per-atom dipole magnitudes
    dm_p = [math.sqrt(sum(d**2 for d in v["dipole"])) for v in mp_p.values()]
    dm_x = [math.sqrt(sum(d**2 for d in v["dipole"])) for v in mp_x.values()]
    # Convert H0 from e·nm to e·Å
    dm_x_A = [m * 10 for m in dm_x]
    print(f"  poltype2 avg per-atom |dipole|: {sum(dm_p)/len(dm_p):.4f} e·Å")
    print(f"  H0 xml   avg per-atom |dipole|: {sum(dm_x_A)/len(dm_x_A):.4f} e·Å")

    # ── QM validation from poltype2 log ──
    print(f"\n[Validation (from poltype2 README)]")
    print(f"  QM dipole  = 13.704 D")
    print(f"  MM dipole  = 13.703 D  (poltype2 AMOEBA, error = 7.3e-5)")
    print(f"  H0 popc.xml: no QM benchmark available (from DOPC prm)")

    print(f"\n[Summary]")
    print(f"  poltype2 final.key : {len(atoms_p)} atom types, QM-derived multipoles,")
    print(f"                       ESP-fitted charges (B3LYP/6-31G*), DB torsions")
    print(f"  H0 popc.xml        : {len(atom_types_x)} atom types, from DOPC.prm,")
    print(f"                       transferable AMOEBA lipid parameters")
    print()


if __name__ == "__main__":
    main()
