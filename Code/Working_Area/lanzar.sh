#!/bin/bash
#SBATCH -t 06:00:00 # execution time. Ex: 6 hours
#SBATCH --mem=5GB
#SBATCH -n 1 -c 4 # number of tasks, number of cores
#SBATCH --mail-type=begin #Envía un correo cuando el trabajo inicia
#SBATCH --mail-type=end #Envía un correo cuando el trabajo finaliza
#SBATCH --mail-user=alejandro.seco@rai.usc.es #Dirección a la que se envía

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

#!/bin/bash
# pip install prolif --quiet
chmod a+rwx Scripts_Carlos/*sh

# Nota: Quitamos -e temporalmente en la llamada del make para controlar
# nosotros mismos el mensaje de error personalizado, o usamos un bloque if.
# Mantenemos las banderas generales de seguridad.
set -uo pipefail

# 1) Muévete al directorio de envío real (Working_Area)
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  cd "$SLURM_SUBMIT_DIR"
fi

DEF_VARS_SCRIPT="Scripts_Carlos/def_vars.sh" #dev_vars.sh or def_vars_mac.sh dependiendo del OS

# =====================================
# Carpeta en Input_Files que va a usarse
# =====================================

INPUT_FOLDER="2REI_Afatinib"
# INPUT_FOLDER="6QH4_Abemaciclib"

# --- Definir variables automáticamente ---
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

# 1. Desactivamos temporalmente el 'set -e' para que el script NO muera si falla make
set +e

# 2. Ejecutamos make all
make -f Makefile_v2 crear_archivos_VR\
    root="$root" \
    output_dir="$output_dir" \
    receptor_pdb="$receptor_pdb" \
    ligand_name="$ligand_name" \
    bfactor="$bfactor" \
    resultados_vina_pdqt="$resultados_vina_pdqt" \
    ff="$ff" \
    water="$water"

# 3. GUARDAMOS inmediatamente si funcionó o falló (0 es éxito, otro número es error)
MAKE_EXIT_CODE=$?

echo "--------------------------------------------------------"
echo "Ejecutando limpieza (make clean)..."

# 3. Hacemos la limpieza, pero si falla, sólo avisamos
make -f Makefile_v2 clean \
    root="$root" \
    output_dir="$output_dir" \
    receptor_pdb="$receptor_pdb" \
    ligand_name="$ligand_name" \
    ff="$ff" \
    water="$water" \
    || echo "Aviso: 'make clean' falló (probablemente no había nada que borrar)."

echo "Limpieza finalizada."
echo "--------------------------------------------------------"

# 4. (Opcional) si quieres volver a tener -e a partir de aquí, lo reactivas
# set -e

# 5. Ahora miramos SOLO el resultado de make all
if [ "$MAKE_EXIT_CODE" -eq 0 ]; then
    echo "Ejecución finalizada correctamente."
    exit 0
else
    echo "!!!! ERROR CRÍTICO !!!!"
    echo "El proceso 'make all' falló (código de error: $MAKE_EXIT_CODE)."
    echo "Se ha realizado el 'make clean', pero revisa los logs de arriba."
    exit 1
fi