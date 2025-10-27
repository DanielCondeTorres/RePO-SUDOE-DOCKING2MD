# RePO-SUDOE-DOCKING2MD

## Install

In CESGA use first:

```
module load cesga/2020 miniconda3/22.11.1-1
```


Use the acpype .yml for Linux

```
conda env create -f nombre_del_archivo.yml
```

## Usage


```
conda activate acype
```

Luego entramos en Working_Area, con cd y damos la ruta a nuestros inputs:


```
 make all root=../experiment_33_docking_files/ output_dir=Dir_name_to_save_results receptor_pdb=receptor.pdb ligand_name=name_of_the_lignand resultados_vina_pdqt=results_out.pdbqt.sdf ff=forcefield_name water=water_model
```
OJOOOO:
ligand_name=Abemaciclib viene de 6TAV_Abemaciclib_out.pdbqt, se tiene que eleminar el 6TAV_ y el _out.pdbqt este make se puede copiar y pegar en el lanzar.sh al final de todo y ejecutarlo con sbatch lanzar.sh

Ex:
```
 make all root=/Users/danielcondetorres/Desktop/GUARDA_PROYECTO/MMDR/Code/Output_Carlos/experiment_33_docking_files/ output_dir=Prueba_final receptor_pdb=3DKO.pdb ligand_name=Abemaciclib  resultados_vina_pdqt=3DKO_Abemaciclib_out.pdbqt.sdf ff=charmm27 water=tip3p
```
## Options
- Water: spc, spce, tip3p, tip4p, tip5p, tips3p
- FF:  charmm36-jul2022, amber99sb-ildn, amber96, charmm27, gromos54a7, gromos53a6, oplsaa

## Necesario:
- Gromacs
- MDanalsys
- mdtraj
- networkx
- prolif


# NOTE!!!!!:

Change in lanzar.sh 

```
make all root=../Input_files/experiment_39_docking_files/ output_dir=Prueba_final receptor_pdb=6QH4.pdb ligand_name=Palbociclib  resultados_vina_pdqt=6QH4_Palbociclib_out.pdbqt.sdf ff=amber99sb-ildn water=tip3p
```

Only is needed to use:

```
sbatch lanzar.sh
```


Necesario darle un par de vueltas al pdbfixer.py, pero es necesario para arreglar los pdb iniciales
