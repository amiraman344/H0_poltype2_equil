# DOPC → POPC AMOEBA Force Field Generation

## Overview

POPC (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine) is a mixed-chain phospholipid with:
- **sn-1 chain**: palmitoyl (C16:0) — fully saturated
- **sn-2 chain**: oleoyl (C18:1, Δ9) — one double bond

No published AMOEBA force field for POPC existed at the time of this work. The approach taken was to derive POPC parameters by combining the published AMOEBA parameters for two closely related lipids — **DOPC** and **POPE** — and creating custom atom types for the sn-1 palmitoyl chain, which differs from both parents.

---

## Source of Published Parameters

| Item | Detail |
|---|---|
| **GitHub repository** | [Inniag/openmm-scripts-amoeba](https://github.com/Inniag/openmm-scripts-amoeba) |
| **Associated paper** | Evaluating Polarizable Biomembrane Simulations against Experiments, *J. Chem. Theory Comput.* **2023**, 19(12), 3956–3970 |
| **DOI** | [10.1021/acs.jctc.3c01333](https://doi.org/10.1021/acs.jctc.3c01333) |
| **PubMed Central** | PMC11137822 |

The repository contains Tinker-format AMOEBA parameter files for several phospholipids. The two files used here are:

| File | Lipid | Tinker FF name |
|---|---|---|
| `dopc.prm` / `dopc.xml` | DOPC (di-oleoyl, C18:1/C18:1) | AMOEBA-DOPC |
| `pope.prm` / `pope.xml` | POPE (palmitoyl-oleoyl-ethanolamine, C16:0/C18:1) | AMOEBA-POPE |

The `.xml` files were converted from Tinker `.prm` format for use with OpenMM. They use a non-standard format (see note on torsions below).

---

## Why DOPC Alone Is Not Enough

DOPC has **two oleoyl chains** (both C18:1). Its sn-1 chain contains a vinyl carbon (C38) adjacent to the double bond, assigned DOPC atom type `1000343` (sp2-adjacent CH₂). In POPC, the sn-1 chain is **palmitoyl** (fully sp3, C16:0), so C38 and all downstream carbons need sp3 CH₂ parameters. Using DOPC types directly causes:

```
Error: Atom C38 of POPC 101 was not assigned to a force field atom type
```

The AMOEBA multipole frame for type `1000343` expects vinyl neighbors in its `kz`/`kx` reference atoms. No such neighbors exist in POPC's palmitoyl chain, so frame assignment fails.

---

## Solution: Custom Types 10004001–10004018

### Strategy

POPE has a palmitoyl sn-1 chain (C16:0) — exactly what POPC needs. The POPE sn-1 parameters use atom types `10002736`–`10002753` (18 types: 9 carbons × 2 atoms each: C and H). These were copied to new custom types `10004001`–`10004018`, with the AMOEBA multipole frames corrected for the different connectivity at the DOPC/POPE junction (C37–C38 bond).

### Type Correspondence

| Custom type | POPE type | Atom(s) in POPC |
|---|---|---|
| 10004001 | 10002736 | C38 (sp3 CH₂) |
| 10004002 | 10002737 | H8X, H8Y |
| 10004003 | 10002738 | C39 |
| 10004004 | 10002739 | H9X, H9Y |
| 10004005 | 10002740 | C310 |
| 10004006 | 10002741 | H10X, H10Y |
| 10004007 | 10002742 | C311 |
| 10004008 | 10002743 | H11X, H11Y |
| 10004009 | 10002744 | C312 |
| 10004010 | 10002745 | H12X, H12Y |
| 10004011 | 10002746 | C313 |
| 10004012 | 10002747 | H13X, H13Y |
| 10004013 | 10002748 | C314 |
| 10004014 | 10002749 | H14X, H14Y |
| 10004015 | 10002750 | C315 |
| 10004016 | 10002751 | H15X, H15Y |
| 10004017 | 10002752 | C316 (methyl terminus) |
| 10004018 | 10002753 | H16X, H16Y, H16Z |

### Frame Reference Fix at C37–C38 Junction

The key correction is in the AMOEBA multipole frame for C38 (type `10004001`). In POPE, C38's frame has:
- `kz` = type `10002734` (POPE C37)
- `kx` = type `10004003` (C39)

In POPC, C37 is a DOPC-type carbon (type `1000341`), not a POPE-type carbon. So:
- `kz` for C38 (custom `10004001`) → `1000341` (DOPC C37)
- `kx` → `10004003` (C39, same custom series)

A second `<Multipole>` entry for the existing DOPC type `1000341` (C37) was also added with `kx="10004001"` so that C37's AMOEBA frame resolves correctly when bonded to the new C38 type.

### Cross-Junction Parameters Added

The C37–C38 bond spans from DOPC parameter space (class `1000323`) to the new custom parameter space (class `10002719`). The following cross-junction terms were added explicitly:

- **1 bond**: `1000323`–`10002719` (C37–C38)
- **4 angles**: C36–C37–C38, C37–C38–C39, H7–C37–C38, C37–C38–H8
- **6 torsions**: all dihedral combinations spanning the C37–C38 bond

---

## Other Fixes Applied to `popc.xml`

### `<Polarize>` tag location

OpenMM 8 requires `<Polarize>` entries to be **inside** `<AmoebaMultipoleForce>`, not in a separate `<AmoebaPolarizeForce>` section. The DOPC-derived XML had them in the wrong location and was corrected.

### Torsion tag name

The published files use `<Angle>` tags *inside* `<AmoebaTorsionForce>` (a non-standard format from the Inniag repository). OpenMM 8.1.1 expects `<Torsion>` tags — the `<Angle>` tags are **silently ignored**. This means torsion terms for the lipid chains are absent from the current `popc.xml`. This is an approximation; proper torsion parameters will come from the poltype2 POPC parameterization (see `POPC_FF_generation/`).

### Malformed v1 torsion entries

Some entries used the old `k1=... k2=... k3=...` attribute style instead of the `amp1`/`angle1` style expected by OpenMM 8. These were removed.

---

## File Descriptions

| File | Description |
|---|---|
| `dopc.prm` | Original Tinker AMOEBA parameter file for DOPC (from Inniag/openmm-scripts-amoeba) |
| `dopc.xml` | OpenMM XML converted from `dopc.prm` |
| `pope.prm` | Original Tinker AMOEBA parameter file for POPE (from Inniag/openmm-scripts-amoeba) |
| `pope.xml` | OpenMM XML converted from `pope.prm` |
| `popc_residue_block.txt` | POPC `<Residue>` XML block with all atom→type assignments |
| `fix_popc_sn1.py` | First-pass script: extracts POPE sn-1 params, creates custom types, patches `popc.xml` |
| `fix_popc_sn1_v2.py` | Revised script with improved section parsing (handles multi-attribute XML tags); reads from `popc_backup.xml` to avoid double-patching |
| `patch_missing.py` | Final patch: adds remaining missing pieces (AtomTypes, VdW, Polarize, torsions) after v2 |
| `popc.xml` | **Final POPC AMOEBA force field** for use with OpenMM (output of the patching pipeline) |

### Script Execution Order

```
fix_popc_sn1_v2.py     ← reads dopc.xml (base), pope.xml (sn-1 source)
                          writes popc.xml (with multipoles, bonds, angles)
        ↓
patch_missing.py        ← reads popc.xml (from above), pope.xml
                          writes final popc.xml (adds AtomTypes, VdW, Polarize, torsions)
```

`fix_popc_sn1.py` is the original first-pass version; `fix_popc_sn1_v2.py` supersedes it.

---

## Limitations and Future Work

1. **Torsion terms absent**: All lipid chain torsion parameters are silently ignored by OpenMM 8.1.1 due to the `<Angle>`-in-`<AmoebaTorsionForce>` format mismatch. This is an approximation.

2. **Poltype2 parameterization**: A proper POPC AMOEBA parameter set is being generated from scratch using poltype2 (Psi4 MP2/6-31G* on POPC fragments). See `../POPC_FF_generation/` for those files. Once complete, `popc.xml` should be replaced with the poltype2-derived `popc_poltype2.xml` for production simulations.

3. **No POPC-specific multipole fitting**: The sn-1 multipoles are transferred from POPE with frame corrections — they are not fitted to POPC's QM electron density. The assumption is that sp3 CH₂ groups in a palmitoyl chain are sufficiently transferable from POPE.
