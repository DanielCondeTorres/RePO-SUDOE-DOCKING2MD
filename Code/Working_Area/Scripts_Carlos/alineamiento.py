#!/usr/bin/env python3
"""
Align a .gro ligand to the first conformer in an .sdf (keeps the original atom order).

Assumptions:
- The SDF is in Ångström (standard). We convert to nm to match GROMACS .gro units.
- The first molecule in the SDF is used (up to the first "$$$$" delimiter).
- The atom ordering between GRO and SDF matches. If atom names match but order doesn't,
  you can use --match-names to reorder GRO by atom name (stable mapping).

Usage:
  python align_gro_to_sdf.py --gro IN.gro --sdf IN.sdf --out OUT.gro [--match-names]

Exit codes:
  0 on success; non-zero if alignment fails.
"""
from __future__ import annotations
import sys
from pathlib import Path
import argparse
from typing import List, Optional

def parse_args(argv: Optional[List[str]]=None):
    p = argparse.ArgumentParser(description="Align a .gro to the first structure in an .sdf (Kabsch).")
    p.add_argument("--gro", required=True, help="Input .gro file (GROMACS units: nm)")
    p.add_argument("--sdf", required=True, help="Input .sdf file (assumed Å). Only first molecule is used.")
    p.add_argument("--out", required=True, help="Output .gro file (aligned)")
    p.add_argument("--match-names", action="store_true",
                   help="If set, reorder GRO atoms to match SDF atom names (stable, first-match).")
    return p.parse_args(argv)

# ---------- Simple parsers ----------

class GRO:
    def __init__(self, title: str, natoms: int, atoms: list, box: str):
        self.title = title
        self.natoms = natoms
        self.atoms = atoms  # list of dicts
        self.box = box

def parse_gro(path: Path) -> GRO:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()

    if len(lines) < 3:
        raise ValueError("GRO file too short.")

    title = lines[0]
    try:
        natoms = int(lines[1].strip())
    except Exception as e:
        raise ValueError(f"Cannot read atom count from GRO: {e}")

    atom_lines = lines[2:2+natoms]
    if len(atom_lines) != natoms:
        raise ValueError("GRO atom lines do not match atom count.")

    atoms = []
    for i, line in enumerate(atom_lines, start=1):
        # GRO fixed-width columns (but many variants exist). Use robust slicing.
        # Typical: resid(5) resname(5) atomname(5) atomnr(5) x(8.3) y(8.3) z(8.3) [vx vy vz]
        resid = line[0:5].strip()
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnr = line[15:20].strip()
        # Positions start around col 20; split rest for robustness
        rest = line[20:].split()
        if len(rest) < 3:
            raise ValueError(f"GRO line {i} missing coordinates: {line}")
        x, y, z = map(float, rest[:3])
        atoms.append({
            "resid": resid,
            "resname": resname,
            "atomname": atomname,
            "atomnr": atomnr,
            "x": x, "y": y, "z": z,
            "raw": line,  # keep original for formatting
        })

    box = lines[2+natoms] if len(lines) >= 3+natoms else ""
    return GRO(title, natoms, atoms, box)

def write_gro(gro: GRO, coords_nm, out_path: Path):
    # coords_nm: list of (x,y,z) in nm, same ordering as gro.atoms
    lines = []
    lines.append(gro.title)
    lines.append(f"{gro.natoms:5d}")
    for atom, (x,y,z) in zip(gro.atoms, coords_nm):
        # Preserve left columns; rewrite coordinates into 8.3f fields
        # Rebuild the line using typical formatting to avoid misalignment
        resid = int(atom["resid"]) if atom["resid"].isdigit() else 0
        resname = atom["resname"]
        atomname = atom["atomname"]
        atomnr = int(atom["atomnr"]) if atom["atomnr"].isdigit() else 0
        lines.append(f"{resid:5d}{resname:>5}{atomname:>5}{atomnr:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    if gro.box.strip():
        lines.append(gro.box)
    else:
        lines.append("   1.00000   1.00000   1.00000")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

class SDF:
    def __init__(self, atoms: List[str], coords_angstrom):
        self.atoms = atoms       # atom symbols or names
        self.coords_A = coords_angstrom  # list of (x,y,z) in Å

def parse_first_sdf_mol(path: Path) -> SDF:
    atoms = []
    coords = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()

    # Find first molecule: counts line is 4th line in V2000
    if len(lines) < 4:
        raise ValueError("SDF file too short / invalid.")
    counts = lines[3]
    if len(counts) < 6:
        raise ValueError("SDF counts line invalid.")
    # Attempt V2000: first 3 ints are natoms, nbonds, ...
    try:
        natoms = int(counts[0:3])
    except ValueError:
        # counts might be space separated
        parts = counts.split()
        natoms = int(parts[0])

    atom_start = 4
    for i in range(natoms):
        line = lines[atom_start + i]
        # V2000 atom block: x y z atomSymbol ...
        try:
            x = float(line[0:10]); y = float(line[10:20]); z = float(line[20:30])
            symbol = line[31:34].strip() or line[30:33].strip() or line.split()[3]
        except Exception:
            parts = line.split()
            x, y, z = map(float, parts[0:3])
            symbol = parts[3] if len(parts) > 3 else "X"
        atoms.append(symbol)
        coords.append((x, y, z))

    # Stop at first $$$$
    # (We already read only the first block via natoms)
    return SDF(atoms, coords)

# ---------- Kabsch alignment ----------

def kabsch(P, Q):
    """
    Compute optimal rotation matrix U and translation t such that:
    Q ~= (U @ P^T).T + t
    where P and Q are Nx3 lists.
    """
    import numpy as np
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    if P.shape != Q.shape or P.shape[1] != 3:
        raise ValueError("P and Q must be Nx3 and same shape")

    # Subtract centroids
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)

    # Covariance
    C = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(C)
    d = (np.linalg.det(V @ Wt) < 0.0)
    if d:
        V[:, -1] *= -1.0
    U = V @ Wt

    # Translation
    t = Q.mean(axis=0) - (U @ P.mean(axis=0))
    return U, t

def apply_transform(P, U, t):
    import numpy as np
    P = np.asarray(P, dtype=float)
    return (U @ P.T).T + t

# ---------- Name matching (optional) ----------

def stable_name_reorder(gro_names: List[str], sdf_names: List[str]) -> List[int]:
    """
    Return indices that map GRO order -> indices in SDF order by atom name (first match, stable).
    If a GRO atom name cannot be found, raise.
    """
    from collections import defaultdict, deque
    pools = defaultdict(deque)
    for i, name in enumerate(gro_names):
        pools[name].append(i)
    order = []
    for name in sdf_names:
        if pools[name]:
            order.append(pools[name].popleft())
        else:
            raise ValueError(f"Could not find GRO atom named '{name}' to match SDF order.")
    return order

# ---------- Main ----------

def main(argv=None):
    args = parse_args(argv)

    gro = parse_gro(Path(args.gro))
    sdf = parse_first_sdf_mol(Path(args.sdf))

    if gro.natoms != len(sdf.coords_A):
        raise SystemExit(f"Atom count mismatch: GRO has {gro.natoms}, SDF has {len(sdf.coords_A)} (first mol).")

    # Prepare coordinate sets
    # GRO coords in nm -> convert to Å for better numeric scale? We'll convert SDF Å -> nm to keep output in nm.
    P_nm = [(a["x"], a["y"], a["z"]) for a in gro.atoms]  # moving (to align)
    Q_nm = [(x*0.1, y*0.1, z*0.1) for (x,y,z) in sdf.coords_A]  # target, in nm

    # Optional: reorder by atom names
    if args.match_names:
        gro_names = [a["atomname"] for a in gro.atoms]
        # SDF has only element symbols; often ACPYPE maps atom names (e.g., C1, H2).
        # If --match-names is used, we try to use element symbols from SDF as names.
        # This only works if your GRO uses element symbols as atom names.
        sdf_names = sdf.atoms
        try:
            idx_map = stable_name_reorder(gro_names, sdf_names)
        except Exception as e:
            raise SystemExit(f"Name matching failed: {e}")
        P_nm = [P_nm[i] for i in idx_map]

    # Compute Kabsch
    U, t = kabsch(P_nm, Q_nm)
    P_aligned_nm = apply_transform(P_nm, U, t)

    # Write output (using original ordering if not reordered; if reordered, we reorder back)
    if args.match_names:
        # Map aligned coords back to original GRO order
        # idx_map: index in SDF order -> index in GRO order? We built it as GRO->SDF, so invert
        inv = [None]*len(P_nm)
        for sdf_idx, gro_idx in enumerate(idx_map):
            inv[gro_idx] = P_aligned_nm[sdf_idx]
        coords_out = inv
    else:
        coords_out = P_aligned_nm

    out_path = Path(args.out)
    write_gro(gro, coords_out, out_path)
    print(f"Wrote aligned GRO: {out_path}")

if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:
        if str(e):
            print(str(e), file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)

