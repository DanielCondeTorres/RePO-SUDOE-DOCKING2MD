#!/bin/bash
#SBATCH -t 06:00:00 # execution time. Ex: 6 hours
#SBATCH --mem-per-cpu=5GB
#SBATCH -n 2 -c 4 # number of tasks, number of cores
##SBATCH --ntasks-per-node=10

# ===============================
# ENTORNO
# ===============================
module load cesga/system miniconda3/22.11.1-1
module load cesga/2020 gcc/system openmpi/4.0.5_ft3_cuda gromacs/2021.4-plumed-2.8.0
conda activate acpype
# pip install prolif --quiet


# ===============================
# CONFIGURACIÓN DEL EXPERIMENTO
# ===============================
chmod a+rwx Scripts_Carlos/*sh
set -euo pipefail

# 1) Muévete al directorio de envío real (Working_Area)
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  cd "$SLURM_SUBMIT_DIR"
fi

# 2) Desde aquí, las rutas relativas son respecto a Working_Area
#    (donde está job.sh y la carpeta Scripts_Carlos)
DEF_VARS_SCRIPT="Scripts_Carlos/def_vars.sh"



# =====================================
# 
# Carpeta en Input_Files que va a usarse
# 
# =====================================

# INPUT_FOLDER="7VON_Alectinib"
INPUT_FOLDER="3P1I_ABEMACICLIB"

# --- Definir variables automáticamente ---
# Aquí se ejecuta tu script definidor y se evalúa su salida
# (asegúrate de que def_vars.sh esté en el mismo directorio que este script)
ROOT="../Input_files/$INPUT_FOLDER"
eval "$("$DEF_VARS_SCRIPT" "$ROOT")"


# Resto de variables fijas
ff="amber99sb-ildn"
water="tip3p"

# ===============================
# EJECUCIÓN
# ===============================
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
echo "Ejecución finalizada correctamente."
