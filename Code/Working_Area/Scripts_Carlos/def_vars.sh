#!/bin/bash
set -euo pipefail

ROOT="${1:-.}"

# --- Buscar receptor: código de 4 caracteres + .pdb (insensible a may/min) ---
RECEPTORS=$(find "$ROOT" -maxdepth 1 -type f \
  -regextype posix-extended \
  -iregex '.*/[A-Za-z0-9]{4}\.pdb' | sort)
num_receptors=$(echo "$RECEPTORS" | grep -c . || true)

if [ "$num_receptors" -eq 0 ]; then
  echo "Error: No se encontró ningún receptor (????.pdb) en '$ROOT'." >&2
  exit 1
elif [ "$num_receptors" -gt 1 ]; then
  echo "Error: Se encontraron múltiples receptores:" >&2
  echo "$RECEPTORS" >&2
  exit 1
fi

receptor_pdb=$(basename "$RECEPTORS")
receptor_code="${receptor_pdb%.pdb}"

# --- Buscar ligando: EXACTAMENTE item1_item2.sdf (solo un guion bajo y sin puntos en item1/2) ---
# Regex POSIX extendido: cada item NO puede contener '/' '_' ni '.'
LIGANDS=$(find "$ROOT" -maxdepth 1 -type f \
  -regextype posix-extended \
  -iregex '.*/[^/_\.]+_[^/_\.]+\.sdf' | sort)
num_ligands=$(echo "$LIGANDS" | grep -c . || true)

if [ "$num_ligands" -eq 0 ]; then
  echo "Error: No se encontró ningún ligando con formato item1_item2.sdf en '$ROOT'." >&2
  exit 1
elif [ "$num_ligands" -gt 1 ]; then
  echo "Error: Se encontraron múltiples ligandos (debe haber solo uno):" >&2
  echo "$LIGANDS" >&2
  exit 1
fi

ligand_base=$(basename "$LIGANDS")
ligand_name="${ligand_base%.sdf}"     # item1_item2 (sin .sdf)
item1="${ligand_name%%_*}"
item2="${ligand_name#*_}"

# Seguridad extra por si algo raro pasa
if [ -z "$item1" ] || [ -z "$item2" ] || echo "$item2" | grep -q '_'; then
  echo "Error: El nombre del ligando no cumple exactamente 'item1_item2.sdf' -> '$ligand_base'." >&2
  exit 1
fi

# --- Construir variables pedidas ---
output_dir="${receptor_code}_${item1}_Out"
resultados_vina_pdqt="${receptor_pdb%.pdb}_${ligand_name}_out.pdbqt.sdf"
bfactor="${receptor_pdb%.pdb}_${ligand_name}_out.pdbqt"

# --- Mostrar en formato listo para 'eval' o 'source' ---
cat <<EOF
root="$ROOT"
output_dir="$output_dir"
receptor_pdb="$receptor_pdb"
ligand_name="$ligand_name"
resultados_vina_pdqt="$resultados_vina_pdqt"
bfactor="$bfactor"
EOF
