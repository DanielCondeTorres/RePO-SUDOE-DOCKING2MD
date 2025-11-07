# import argparse
# import prolif as plf
# import MDAnalysis as mda
# import matplotlib.pyplot as plt
# from prolif.plotting.network import LigNetwork

# # ------------------ PARSER DE ARGUMENTOS ------------------
# parser = argparse.ArgumentParser(description="Genera red de interacciones proteína-ligando con ProLIF")
# parser.add_argument("--protein", required=True, help="Ruta al archivo PDB de la proteína")
# parser.add_argument("--poses", required=True, help="Ruta al archivo SDF con las poses de ligandos")
# parser.add_argument("--output", default="interacciones.html", help="Ruta de salida para la imagen (PNG/SVG/PDF)")
# parser.add_argument("--pose_index", type=int, default=0, help="Índice de la pose del ligando a analizar (default=0)")
# parser.add_argument("--residue_shift", type=int, default=0, help="Valor constante a sumar a los IDs de los residuos")
# args = parser.parse_args()

# # ------------------ CARGAR Y MODIFICAR PROTEÍNA ------------------
# protein = mda.Universe(args.protein)
# protein_atoms = protein.select_atoms("protein")

# # Ajustar IDs de residuos
# for residuo in protein_atoms.residues:
#     residuo.resid += args.residue_shift

# # Crear objeto ProLIF para la proteína
# protein_plf = plf.Molecule.from_mda(protein,NoImplicit=False)

# # ------------------ CARGAR LIGANDOS ------------------
# poses_plf = plf.sdf_supplier(args.poses)

# # ------------------ CALCULAR INTERACCIONES ------------------
# fp = plf.Fingerprint(count=True, vicinity_cutoff=6.0)
# fp.run_from_iterable(poses_plf, protein_plf)

# # ------------------ CREAR Y GUARDAR GRÁFICO ------------------
# lignetwork = LigNetwork.from_fingerprint(
#     fp,
#     ligand_mol=poses_plf[args.pose_index],
#     kind="frame",  
#     frame=args.pose_index,
#     display_all=False,
#     threshold=0.6
# )

# lignetwork.save(args.output)

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
print("🔄 Cargando proteína...")
protein = mda.Universe(args.protein)

# *[FIX APLICADO]* Seleccionar solo los átomos de la proteína ANTES de adivinar bonds
protein_atoms = protein.select_atoms("protein")
try:
    # Intenta adivinar bonds solo en los átomos de la proteína
    protein_atoms.guess_bonds()  
    print("✅ Bonds detectados correctamente en la proteína.")
except ValueError as e:
    print(f"⚠️ Advertencia: Error al adivinar bonds. Intenta un método alternativo.")
    print(f"Detalle del error: {e}")
    # Si falla, se puede intentar adivinar bonds en el universo completo con un vdw radius para 'CO'
    if 'vdw radii for types: CO' in str(e):
        print("💡 Aplicando workaround para el tipo 'CO' con vdw_radii={'CO': 1.8}.")
        vdw_radii = {'CO': 1.8}
        protein.atoms.guess_bonds(vdwradii=vdw_radii)
    else:
        # Si es otro error, lanzarlo
        raise e

# Ajustar IDs de residuos
for residuo in protein_atoms.residues:
    residuo.resid += args.residue_shift
# Crear objeto ProLIF para la proteína
protein_plf = plf.Molecule.from_mda(protein, NoImplicit=False)

# ------------------ CARGAR LIGANDOS ------------------
print("🔄 Cargando poses...")
poses_plf = plf.sdf_supplier(args.poses)

# ------------------ CALCULAR INTERACCIONES ------------------
print("🔄 Calculando interacciones...")
fp = plf.Fingerprint(count=True, vicinity_cutoff=6.0)
fp.run_from_iterable(poses_plf, protein_plf)

# ------------------ DEBUG: Verificar interacciones ------------------
df = fp.to_dataframe()
print(f"📊 Interacciones detectadas: {df.shape[0]} filas")
print(f"   Columnas: {df.columns.tolist()}")
if not df.empty:
    print("   Ejemplo:\n", df.head(3))
else:
    print("❌ ¡NO HAY INTERACCIONES!")
    print("💡 Prueba:")
    print("   - --pose_index X (X=mejor pose del docking)")
    print("   - vicinity_cutoff=10.0 en Fingerprint")
    print("   - Verifica distancia ligando-proteína en PyMOL")
    exit(1)

# ------------------ CREAR Y GUARDAR GRÁFICO ------------------
print("🎨 Generando red...")
lignetwork = LigNetwork.from_fingerprint(
    fp,
    ligand_mol=poses_plf[args.pose_index],
    kind="frame",
    frame=args.pose_index,
    display_all=False,
    threshold=0.6
)
lignetwork.save(args.output)
print(f"✅ ¡Red guardada en: {args.output}!")
print("🔗 Ábrela en tu navegador para ver las interacciones interactivas.")