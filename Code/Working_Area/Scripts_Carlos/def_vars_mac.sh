#!/bin/bash
set -euo pipefail

ROOT="${1:-.}"

# Inicializa variables para evitar “unbound variable”
RECEPTORS=""
LIGANDS=""
receptor_pdb=""
receptor_code=""
ligand_base=""
ligand_name=""
item1=""
item2=""

# --- Buscar receptor ----
RECEPTORS=$(find "$ROOT" -maxdepth 1 -type f -iname '*.pdb' \
  -exec basename {} \; | grep -E '^[A-Za-z0-9]{4}\.pdb$' | sort || true)
num_receptors=$(echo "$RECEPTORS" | grep -c . || true)
if [ "$num_receptors" -eq 0 ]; then
  echo "Error: no receptor (????.pdb) en $ROOT" >&2
  exit 1
elif [ "$num_receptors" -gt 1 ]; then
  echo "Error: múltiples receptores:" >&2
  echo "$RECEPTORS" >&2
  exit 1
fi
receptor_pdb="$RECEPTORS"
receptor_code="${receptor_pdb%.pdb}"

# --- Buscar ligando ----
LIGANDS=$(find "$ROOT" -maxdepth 1 -type f -iname '*.sdf' \
  -exec basename {} \; | grep -E '^[^/_\.]+_[^/_\.]+\.sdf$' | sort || true)
num_ligands=$(echo "$LIGANDS" | grep -c . || true)
if [ "$num_ligands" -eq 0 ]; then
  echo "Error: no ligando item1_item2.sdf en $ROOT" >&2
  exit 1
elif [ "$num_ligands" -gt 1 ]; then
  echo "Error: múltiples ligandos:" >&2
  echo "$LIGANDS" >&2
  exit 1
fi
ligand_base=$(basename "$LIGANDS")
ligand_name="${ligand_base%.sdf}"
item1="${ligand_name%%_*}"
item2="${ligand_name#*_}"
if [ -z "$item1" ] || [ -z "$item2" ] || echo "$item2" | grep -q '_'; then
  echo "Error: nombre de ligando no cumple item1_item2.sdf -> $ligand_base" >&2
  exit 1
fi

# --- Construir variables ---
output_dir="${receptor_code}_${item1}_Out"
resultados_vina_pdqt="${receptor_code}_${ligand_name}_out.pdbqt.sdf"
bfactor="${receptor_code}_${ligand_name}_out.pdbqt"

# --- Imprimir en formato para source ---
cat <<EOF
root="$ROOT"
output_dir="$output_dir"
receptor_pdb="$receptor_pdb"
ligand_name="$ligand_name"
bfactor="$bfactor"
resultados_vina_pdqt="$resultados_vina_pdqt"
EOF

