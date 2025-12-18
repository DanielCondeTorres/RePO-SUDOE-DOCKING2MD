#!/bin/bash
module load cesga/2020 gcc/system openmpi/4.0.5_ft3_cuda gromacs/2021.4-plumed-2.8.0

# Número de cores disponibles para esta tarea (lo da SLURM). Si no está definido, asumimos 4.
NCORES=${SLURM_CPUS_PER_TASK:-4}

# Función envoltorio para lanzar mdrun intentando varias combinaciones de ntmpi/ntomp.
# Usa una corrida corta con -nsteps 0 sólo para probar la domain decomposition.
run_mdrun_safe() {
    # $1 = deffnm (nombre base de los archivos de salida)
    # $@ (a partir de $2) = argumentos extra para mdrun (-v, -cpi, -noappend, etc.)
    local DEFFNM="$1"; shift
    local EXTRA_ARGS=("$@")

    # Filtramos opciones que NO queremos usar en la corrida de prueba
    # (por ejemplo -cpi y -noappend, que pueden fallar si no hay checkpoint todavía).
    local TEST_ARGS=()
    for arg in "${EXTRA_ARGS[@]}"; do
        case "$arg" in
            -cpi|-noappend)
                ;; # las omitimos en el test
            *)
                TEST_ARGS+=("$arg")
                ;;
        esac
    done

    # Candidatos de ntmpi: primero más agresivo, luego más conservador
    local CANDIDATES=(4 2 1)

    for NMPI in "${CANDIDATES[@]}"; do
        # No usar más ranks de los cores disponibles
        if (( NMPI > NCORES )); then
            continue
        fi

        # Hilos OpenMP por rank (como mínimo 1)
        local NTOMP=$(( NCORES / NMPI ))
        if (( NTOMP < 1 )); then
            NTOMP=1
        fi

        echo ">>> Probando configuración DD: ntmpi=${NMPI}, ntomp=${NTOMP} (deffnm=${DEFFNM}, test con -nsteps 0)"

        # Corrida de prueba: sólo setup + DD, 0 pasos de integración
        gmx mdrun -deffnm "${DEFFNM}" -ntmpi "${NMPI}" -ntomp "${NTOMP}" -nsteps 0 \
                  "${TEST_ARGS[@]}" &> "${DEFFNM}_ddtest.log"

        if [ $? -eq 0 ]; then
            echo ">>> OK con ntmpi=${NMPI}, ntomp=${NTOMP}. Lanzando producción real..."
            export OMP_NUM_THREADS="${NTOMP}"
            gmx mdrun -deffnm "${DEFFNM}" -ntmpi "${NMPI}" -ntomp "${NTOMP}" \
                      "${EXTRA_ARGS[@]}"
            return $?
        else
            echo ">>> FALLÓ DD con ntmpi=${NMPI}, ntomp=${NTOMP}. Probando siguiente..."
        fi
    done

    echo "ERROR: No se encontró ninguna configuración válida de ntmpi/ntomp para ${DEFFNM}" >&2
    return 1
}

COMPLEX_PDB_FILE_1="$1"  # Archivo PDB del complejo
OUTPUT_SAVE="$2"

gmx solvate -cp "$COMPLEX_PDB_FILE_1" -cs spc216.gro -o "$OUTPUT_SAVE"complex_water_box.pdb -p "$OUTPUT_SAVE"topol.top

# Add ions
gmx grompp -f ../Input_files/MDP_FILES/ions.mdp -c "$OUTPUT_SAVE"complex_water_box.pdb -p "$OUTPUT_SAVE"topol.top -o "$OUTPUT_SAVE"ions.tpr -maxwarn 1
echo "SOL" | gmx genion -s "$OUTPUT_SAVE"ions.tpr -o "$OUTPUT_SAVE"complex_solv_ions.gro -p "$OUTPUT_SAVE"topol.top -pname NA -nname CL -neutral

# Minimization of energy
gmx grompp -f ../Input_files/MDP_FILES/min.mdp -c "$OUTPUT_SAVE"complex_solv_ions.gro -p "$OUTPUT_SAVE"topol.top -o "$OUTPUT_SAVE"min_fep1.tpr -maxwarn 2
run_mdrun_safe "${OUTPUT_SAVE}min_fep1" -v || exit 1

# Equilibration NVT
gmx grompp -f ../Input_files/MDP_FILES/eq_nvt.mdp -c "$OUTPUT_SAVE"min_fep1.gro -r "$OUTPUT_SAVE"min_fep1.gro -p "$OUTPUT_SAVE"topol.top -o "$OUTPUT_SAVE"eq_nvt_fep.tpr -maxwarn 2
export GMX_MAXCONSTRWARN=-1
run_mdrun_safe "${OUTPUT_SAVE}eq_nvt_fep" -v || exit 1
unset GMX_MAXCONSTRWARN

# Equilibración NPT 1
gmx grompp -f ../Input_files/MDP_FILES/eq_npt_fep_1.mdp -c "$OUTPUT_SAVE"eq_nvt_fep.gro -p "$OUTPUT_SAVE"topol.top -o "$OUTPUT_SAVE"eq_npt_fep_1.tpr -r "$OUTPUT_SAVE"eq_nvt_fep.gro -maxwarn 2
run_mdrun_safe "${OUTPUT_SAVE}eq_npt_fep_1" -v || exit 1

#Create production file: prod.tpr
gmx grompp -f ../Input_files/MDP_FILES/prod.mdp -c "$OUTPUT_SAVE"eq_npt_fep_1.gro -p "$OUTPUT_SAVE"topol.top -o "$OUTPUT_SAVE"prod.tpr -maxwarn 2

# Si se hace en gromacs igual es mejor lanzar un .sh
#srun gmx_mpi mdrun -pin on -cpi -noappend -s "$OUTPUT_SAVE"prod.tpr
run_mdrun_safe "${OUTPUT_SAVE}prod" -v -cpi -noappend || exit 1


#Añado analisis aqui:
bash Scripts_Carlos/analyisis_creation.sh "$OUTPUT_SAVE"