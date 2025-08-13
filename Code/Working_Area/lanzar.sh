#!/bin/bash
<<<<<<< HEAD
#SBATCH -t  6:00:00 # execution time. Ex: 1 hour
=======
#SBATCH -t  2:00:00 # execution time. Ex: 1 hour
>>>>>>> 887bb383f929e7141fac86b4df57eac87cef17c3
#SBATCH --mem-per-cpu=1GB
#SBATCH -n 10 -c 5# number of tasks, number of cores
#SBATCH --ntasks-per-node=10
module load cesga/system miniconda3/22.11.1-1
module load cesga/2020 gcc/system openmpi/4.0.5_ft3_cuda gromacs/2021.4-plumed-2.8.0
conda activate acpype
make all root=../Input_files/experiment_39_docking_files/ output_dir=Prueba_final receptor_pdb=6QH4.pdb ligand_name=Palbociclib  resultados_vina_pdqt=6QH4_Palbociclib_out.pdbqt.sdf ff=amber99sb-ildn water=tip3p
