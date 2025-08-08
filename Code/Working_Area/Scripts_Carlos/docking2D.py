import argparse
import prolif as plf
import MDAnalysis as mda
import matplotlib.pyplot as plt
from prolif.plotting.network import LigNetwork

# ------------------ PARSER DE ARGUMENTOS ------------------
parser = argparse.ArgumentParser(description="Genera red de interacciones proteína-ligando con ProLIF")
parser.add_argument("--protein", required=True, help="Ruta al archivo PDB de la proteína")
parser.add_argument("--poses", required=True, help="Ruta al archivo SDF con las poses de ligandos")
parser.add_argument("--output", default="interacciones.html", help="Ruta de salida para la imagen (PNG/SVG/PDF)")
parser.add_argument("--pose_index", type=int, default=0, help="Índice de la pose del ligando a analizar (default=0)")
parser.add_argument("--residue_shift", type=int, default=0, help="Valor constante a sumar a los IDs de los residuos")
args = parser.parse_args()

# ------------------ CARGAR Y MODIFICAR PROTEÍNA ------------------
protein = mda.Universe(args.protein)
protein_atoms = protein.select_atoms("protein")

# Ajustar IDs de residuos
for residuo in protein_atoms.residues:
    residuo.resid += args.residue_shift

# Crear objeto ProLIF para la proteína
protein_plf = plf.Molecule.from_mda(protein,NoImplicit=False)

# ------------------ CARGAR LIGANDOS ------------------
poses_plf = plf.sdf_supplier(args.poses)

# ------------------ CALCULAR INTERACCIONES ------------------
fp = plf.Fingerprint(count=True, vicinity_cutoff=6.0)
fp.run_from_iterable(poses_plf, protein_plf)

# ------------------ CREAR Y GUARDAR GRÁFICO ------------------
lignetwork = LigNetwork.from_fingerprint(
    fp,
    ligand_mol=poses_plf[args.pose_index],
    kind="frame",  
    frame=args.pose_index,
    display_all=False,
    threshold=0.6
)

lignetwork.save(args.output)

