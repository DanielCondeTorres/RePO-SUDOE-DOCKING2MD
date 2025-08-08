import prolif as plf
import MDAnalysis as mda
import json
import re

# --- Configuración inicial ---
file     = '../../Output'
pdb_path = f"{file}/receptor.pdb"
sdf_path = f"{file}/all.sdf"

# Carga de la proteína y ajuste de numeración (opcional)
protein = mda.Universe(pdb_path)
for residuo in protein.select_atoms("protein").residues:
    residuo.resid += 0  # ajusta si necesitas un offset

# Creamos el objeto ProLIF
protein_plf = plf.Molecule.from_mda(protein)
poses_plf    = plf.sdf_supplier(sdf_path)

# Inicializar y calcular el fingerprint
fp = plf.Fingerprint(count=True, vicinity_cutoff=6.0)
fp.run_from_iterable(poses_plf, protein_plf)

# Regex para parsear "ALA1148.A"
pattern = re.compile(r"^([A-Z]{3})(\d+)\.([A-Za-z0-9])$")

# --------------------------------------------------
# 1) Leer todas las energías VINA RESULT del SDF
# --------------------------------------------------
vina_re = re.compile(r"VINA RESULT:\s*([-+]?\d*\.\d+|\d+)")
energies = []
with open(sdf_path, 'r') as f:
    for line in f:
        m = vina_re.search(line)
        if m:
            energies.append(float(m.group(1)))

# --------------------------------------------------
# 2) Preparamos el DataFrame completo de interacciones
# --------------------------------------------------
df_full = fp.to_dataframe().reset_index()

# --------------------------------------------------
# 3) Iteramos sobre cada pose para armar el JSON
# --------------------------------------------------
results = []
N_poses = len(poses_plf)
for pose_index in range(N_poses):
    # Extraemos la energía (si existe)
    energy = energies[pose_index] if pose_index < len(energies) else None

    # Filtramos interacciones de este frame
    df_pose = df_full[df_full["Frame"] == pose_index]
    row     = df_pose.iloc[0]

    # Identificamos residuos con valor > 0
    interacting = {
        col[1]
        for col, val in row.items()
        if isinstance(col, tuple) and val > 0
    }

    # Parseamos cada string de residuo a dict base
    residues = []
    for entry in sorted(interacting):
        m = pattern.match(entry)
        if not m:
            continue
        resname, resid, chain = m.groups()
        residues.append({
            "resname": resname,
            "resid": int(resid),
            "chain": chain
        })

    # Para cada residuo, añadimos la lista de índices de sus átomos
    for res in residues:
        sel = (f"resname {res['resname']} "
               f"and resid {res['resid']} "
               f"and segid {res['chain']}")
        atom_indices = protein.select_atoms(sel).indices.tolist()
        res["atoms"] = atom_indices

    # Guardamos el bloque de esta pose
    results.append({
        "model": float(pose_index)+1,
        "energy": energy,
        "residues": residues
    })
import json

# Assuming `results` is the list of dictionaries from the previous script
# and `file` is the directory path.
file_directory = '../../Output'
output_path = f"{file_directory}/aa_interactions_arx.json"

# Example `results` placeholder (remove if running with real `results`)
# results = [...]

# Save JSON to the specified file
with open(output_path, 'w') as json_file:
    json.dump(results, json_file, indent=2)

print(f"JSON de interacciones guardado en: {output_path}")
# --------------------------------------------------
# 4) Imprimimos el JSON final
# --------------------------------------------------
print(json.dumps(results, indent=2))
