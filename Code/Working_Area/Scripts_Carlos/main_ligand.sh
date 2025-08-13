#!/bin/bash
module load cesga/2020 gcc/system openmpi/4.0.5_ft3_cuda gromacs/2021.4-plumed-2.8.0

COMPLEX_PDB_FILE_1="$1"  # Archivo PDB del complejo
OUTPUT_SAVE="$2"
gmx solvate -cp "$COMPLEX_PDB_FILE_1" -cs spc216.gro -o "$OUTPUT_SAVE"complex_water_box.pdb -p "$OUTPUT_SAVE"topol.top
# Add ions


gmx grompp -f ../Input_files/MDP_FILES/ions.mdp -c "$OUTPUT_SAVE"complex_water_box.pdb -p "$OUTPUT_SAVE"topol.top -o "$OUTPUT_SAVE"ions.tpr -maxwarn 1
echo "SOL" | gmx genion -s "$OUTPUT_SAVE"ions.tpr -o "$OUTPUT_SAVE"complex_solv_ions.gro -p "$OUTPUT_SAVE"topol.top -pname NA -nname CL -neutral
# Minimization of energy
gmx grompp -f ../Input_files/MDP_FILES/min.mdp -c "$OUTPUT_SAVE"complex_solv_ions.gro -p "$OUTPUT_SAVE"topol.top -o "$OUTPUT_SAVE"min_fep1.tpr -maxwarn 2
gmx mdrun -v -deffnm "$OUTPUT_SAVE"min_fep1 -ntmpi 4

# Equilibration NVT
gmx grompp -f ../Input_files/MDP_FILES/eq_nvt.mdp -c "$OUTPUT_SAVE"min_fep1.gro -r "$OUTPUT_SAVE"min_fep1.gro -p "$OUTPUT_SAVE"topol.top -o "$OUTPUT_SAVE"eq_nvt_fep.tpr -maxwarn 2
gmx mdrun -v -deffnm "$OUTPUT_SAVE"eq_nvt_fep -ntmpi 4
# Equilibración NPT 1
gmx grompp -f ../Input_files/MDP_FILES/eq_npt_fep_1.mdp -c "$OUTPUT_SAVE"eq_nvt_fep.gro -p "$OUTPUT_SAVE"topol.top -o "$OUTPUT_SAVE"eq_npt_fep_1.tpr -r "$OUTPUT_SAVE"eq_nvt_fep.gro -maxwarn 2
gmx mdrun -v -deffnm "$OUTPUT_SAVE"eq_npt_fep_1 -ntmpi 4
#Create production file: prod.tpr
gmx grompp -f ../Input_files/MDP_FILES/prod.mdp -c "$OUTPUT_SAVE"eq_npt_fep_1.gro -p "$OUTPUT_SAVE"topol.top -o "$OUTPUT_SAVE"prod.tpr -maxwarn 2

# Si se hace en gromacs igual es mejor lanzar un .sh
#srun gmx_mpi mdrun -pin on -cpi -noappend -s "$OUTPUT_SAVE"prod.tpr
gmx mdrun -v -deffnm "$OUTPUT_SAVE"prod -ntmpi 4


#Añado analisis aqui:
bash Scripts_Carlos/analyisis_creation.sh "$OUTPUT_SAVE"
