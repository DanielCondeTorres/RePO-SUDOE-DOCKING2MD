# RePO-SUDOE-DOCKING2MD

## Usage

```
 make all root=../experiment_33_docking_files/ output_dir=Dir_name_to_save_results receptor_pdb=receptor.pdb ligand_name=name_of_the_lignand resultados_vina_pdqt=results_out.pdbqt.sdf ff=forcefield_name water=water_model
```

Ex:
```
make all root=../experiment_33_docking_files/ output_dir=Prueba receptor_pdb=3DKO.pdb ligand_name=Abemaciclib  resultados_vina_pdqt=3DKO_Abemaciclib_out.pdbqt.sdf ff=amber99sb-ildn water=tip3p

```
## Options
- Water: spc, spce, tip3p, tip4p, tip5p, tips3p
- FF:  charmm36-jul2022, amber99sb-ildn, amber96, charmm27, gromos54a7, gromos53a6, oplsaa

## Necesario:
- Gromacs
- MDanalsys
- mdtraj
