#!/usr/bin/env python3
"""
Align a moving PDB onto a reference PDB using protein atoms (CA by default)
with an ICP (Iterative Closest Point) approach + Kabsch rotation.

It finds the optimal rigid-body transform even if the proteins have
different atom counts or residue mismatches, then applies that transform
to ALL atoms in the moving PDB — preserving connectivity.

Usage:
  python align_pdb_protein_icp.py \
      --ref 2REI.pdb --mov receptor_clean.pdb --out receptor_aligned.pdb

Options:
  --backbone-only   Use N, CA, C, O atoms for alignment instead of CA-only
  --max-dist 4.0    Maximum distance (Å) for matching during ICP
  --max-iter 80     Maximum ICP iterations
  --tol 1e-4        RMSD convergence tolerance
"""

from __future__ import annotations
import sys
from pathlib import Path
import argparse
import numpy as np
from typing import List, Tuple

# ---------------- CLI ----------------

def parse_args():
    p = argparse.ArgumentParser(description="Align two PDBs using ICP (protein-only rigid-body alignment).")
    p.add_argument("--ref", required=True, help="Reference PDB (target).")
    p.add_argument("--mov", required=True, help="Moving PDB (will be aligned).")
    p.add_argument("--out", required=True, help="Output aligned PDB.")
    p.add_argument("--backbone-only", action="store_true", help="Use backbone atoms (N, CA, C, O) instead of CA-only.")
    p.add_argument("--max-dist", type=float, default=4.0, help="Max distance (Å) for matching (default 4.0).")
    p.add_argument("--max-iter", type=int, default=80, help="Max ICP iterations (default 80).")
    p.add_argument("--tol", type=float, default=1e-4, help="Convergence tolerance (default 1e-4).")
    return p.parse_args()

# ---------------- PDB parsing ----------------

AMINO = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "HSD","HSE","HSP","HIE","HID","HIP","MSE","SEC","PYL",
    "PCA","CSD","CSO","CME","CYS2","ACE","NME"
}
WATER = {"HOH","WAT","H2O","SOL","TIP3","TIP4","TIP5","SPC","SPCE"}
BACKBONE = {"N","CA","C","O"}

class Atom:
    __slots__ = ("i","record","name","resname","chain","resseq","altloc","x","y","z","line")
    def __init__(self,i,record,name,resname,chain,resseq,altloc,x,y,z,line):
        self.i=i; self.record=record; self.name=name; self.resname=resname
        self.chain=chain; self.resseq=resseq; self.altloc=altloc
        self.x=x; self.y=y; self.z=z; self.line=line

def parse_pdb(path: Path):
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms=[]
    for idx,line in enumerate(lines):
        rec=line[0:6]
        if rec not in ("ATOM  ","HETATM"):
            continue
        name=line[12:16].strip()
        altloc=line[16:17].strip()
        resname=line[17:20].strip()
        chain=line[21:22].strip()
        resseq=line[22:26].strip()
        try:
            x=float(line[30:38]); y=float(line[38:46]); z=float(line[46:54])
        except Exception:
            parts=line.split()
            if len(parts)>=3:
                x,y,z=map(float,parts[-3:])
            else:
                continue
        atoms.append(Atom(idx,rec,name,resname,chain,resseq,altloc,x,y,z,line))
    return lines, atoms

def is_protein_like(a: Atom):
    if a.resname in WATER: return False
    if a.record=="ATOM  ": return True
    if a.record=="HETATM" and a.resname in AMINO: return True
    return False

def select_protein_atoms(atoms: List[Atom], backbone_only=False):
    sel = [a for a in atoms if is_protein_like(a)]
    if backbone_only:
        sel = [a for a in sel if a.name in BACKBONE]
    else:
        sel = [a for a in sel if a.name == "CA"]
    return sel

# ---------------- Math: Kabsch + ICP ----------------

def kabsch(P,Q):
    P=np.asarray(P,float); Q=np.asarray(Q,float)
    Pc=P-P.mean(axis=0); Qc=Q-Q.mean(axis=0)
    C=Pc.T@Qc
    V,S,Wt=np.linalg.svd(C)
    if np.linalg.det(V@Wt)<0: V[:,-1]*=-1
    U=V@Wt
    t=Q.mean(axis=0) - U@P.mean(axis=0)
    return U,t

def apply_transform(P,U,t):
    return (U@P.T).T + t

def pair_greedy_unique(P,Q,max_dist):
    used=np.zeros(len(Q),dtype=bool)
    Ii=[]; Jj=[]
    for i in range(len(P)):
        d2=((Q-P[i])**2).sum(axis=1)
        order=np.argsort(d2)
        for j in order:
            if used[j]: continue
            if d2[j] <= max_dist**2:
                used[j]=True
                Ii.append(i); Jj.append(j)
                break
            else:
                break
    return Ii,Jj

def rmsd(A,B):
    if len(A)==0: return float("inf")
    A=np.asarray(A,float); B=np.asarray(B,float)
    return np.sqrt(((A-B)**2).sum()/len(A))

def icp(P0,Q,max_dist=4.0,max_iter=80,tol=1e-4):
    P=np.array(P0,float)
    prev=float("inf")
    for it in range(1,max_iter+1):
        Ii,Jj=pair_greedy_unique(P,Q,max_dist)
        if len(Ii)<3: break
        Pm,Pn=P[Ii],Q[Jj]
        U,t=kabsch(Pm,Pn)
        P=apply_transform(P,U,t)
        r=rmsd(Pm,Pn)
        print(f"Iter {it:02d}: matches={len(Ii)} RMSD={r:.3f} Å")
        if abs(prev-r)<tol: break
        prev=r
    U_final,t_final=kabsch(P0,P)
    return U_final,t_final,P

# ---------------- Writer ----------------

def write_pdb(lines: List[str], atoms: List[Atom], U, t, out_path: Path):
    new_lines=list(lines)
    for atom in atoms:
        x2,y2,z2=(U@np.array([atom.x,atom.y,atom.z]))+t
        line=new_lines[atom.i]
        if len(line)<54: line=line.ljust(54)
        xyz=f"{x2:8.3f}{y2:8.3f}{z2:8.3f}"
        new_lines[atom.i]=line[:30]+xyz+line[54:]
    content="\n".join(new_lines)
    if not content.endswith("\n"): content+="\n"
    out_path.write_text(content,encoding="utf-8")

# ---------------- Main ----------------

def main():
    args=parse_args()
    ref_lines,ref_atoms=parse_pdb(Path(args.ref))
    mov_lines,mov_atoms=parse_pdb(Path(args.mov))

    ref_sel=select_protein_atoms(ref_atoms, backbone_only=args.backbone_only)
    mov_sel=select_protein_atoms(mov_atoms, backbone_only=args.backbone_only)

    Q=np.array([(a.x,a.y,a.z) for a in ref_sel])
    P0=np.array([(a.x,a.y,a.z) for a in mov_sel])

    print(f"Selected {len(ref_sel)} atoms in REF, {len(mov_sel)} atoms in MOV for alignment.")
    U,t,P_aligned = icp(P0,Q,max_dist=args.max_dist,max_iter=args.max_iter,tol=args.tol)

    write_pdb(mov_lines,mov_atoms,U,t,Path(args.out))
    print(f"\nAligned PDB written to: {args.out}")

if __name__=="__main__":
    try:
        main()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
