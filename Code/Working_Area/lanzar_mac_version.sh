#!/bin/bash
# ==========================================
#   Script principal de ejecución
# ==========================================
HISTTIMEFORMAT=${HISTTIMEFORMAT-""}
set -euo pipefail

chmod a+rwx Scripts_Carlos/*sh

# Si se lanza desde SLURM, moverse al directorio correcto
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  cd "$SLURM_SUBMIT_DIR"
fi

DEF_VARS_SCRIPT="Scripts_Carlos/def_vars_mac.sh"
INPUT_FOLDER="experiment_105_test"
ROOT="../Input_files/$INPUT_FOLDER"

# -------------------------------
# Ejecutar script de variables
# -------------------------------
if ! vars_output="$("$DEF_VARS_SCRIPT" "$ROOT")"; then
  echo "Error: fallo al ejecutar $DEF_VARS_SCRIPT" >&2
  exit 1
fi

# Cargar variables definidas (sin eval)
tempfile=$(mktemp)
echo "$vars_output" > "$tempfile"
# shellcheck disable=SC1090
source "$tempfile"
rm -f "$tempfile"

# -------------------------------
# Parámetros fijos
# -------------------------------
ff="amber99sb-ildn"
water="tip3p"

# -------------------------------
# Mostrar resumen
# -------------------------------
echo "--------------------------------------------------------"
echo "Iniciando experimento con los siguientes parámetros:"
echo "root                = $root"
echo "output_dir          = $output_dir"
echo "receptor_pdb        = $receptor_pdb"
echo "ligand_name         = $ligand_name"
echo "bfactor             = $bfactor"
echo "resultados_vina_pdqt= $resultados_vina_pdqt"
echo "ff                  = $ff"
echo "water               = $water"
echo "--------------------------------------------------------"

# -------------------------------
# Ejecutar make
# -------------------------------
make all \
  root="$root" \
  output_dir="$output_dir" \
  receptor_pdb="$receptor_pdb" \
  ligand_name="$ligand_name" \
  bfactor="$bfactor" \
  resultados_vina_pdqt="$resultados_vina_pdqt" \
  ff="$ff" \
  water="$water"

echo "--------------------------------------------------------"
echo "Ejecución completada correctamente."
echo "--------------------------------------------------------"

