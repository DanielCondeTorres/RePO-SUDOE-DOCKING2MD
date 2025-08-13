SCRIPT = main.sh



# make all root=/Users/danielcondetorres/Desktop/GUARDA_PROYECTO/MMDR/Code/Output_Carlos/experiment_33_docking_files/ output_dir=Prueba receptor_pdb=3DKO.pdb ligand_name=Abemaciclib  resultados_vina_pdqt=3DKO_Abemaciclib_out.pdbqt.sdf ff=amber99sb-ildn water=tip3p
.PHONY: all steps_previos_al_docking_proteina steps_previos_al_docking_ligando docking2D_resutls combinar_topologias obtener_pdb_combinado agua_e_iones analysis crear_archivos_VR clean
all: steps_previos_al_docking_proteina steps_previos_al_docking_ligando docking2D_resutls combinar_topologias obtener_pdb_combinado agua_e_iones analysis crear_archivos_VR clean


# Pedir a Carlos que use estos pdb para hacer el docking
# output_dir: carpeta donde se va a guardar la simulacion
# root: carpeta donde se encuentran ACTUALMENTE los archivos de la web de Carlos .sdf, .pdb ... Los resultados del Docking vamos
# receptor_pdb: .pdb del receptor a emplear
# ligand_name: Nombre del ligando a emplear (porque los .pdb, .sdf, y la carpeta que se va a crear se llaman todas igual)

#make steps_previos_al_docking_proteina root=/Users/danielcondetorres/Desktop/GUARDA_PROYECTO/MMDR/Code/Output_Carlos/experiment_33_docking_files/ output_dir=Prueba receptor_pdb=3DKO.pdb ff=amber99sb-ildn water=tip3p
steps_previos_al_docking_proteina:
	mkdir -p  ../$(output_dir)
	conda run -n acpype python Scripts_Carlos/fixer.py -i $(root)/$(receptor_pdb) -o ../$(output_dir)/receptor_clean_previo.pdb
	cp $(root)/$(receptor_pdb) ../$(output_dir)/
	gmx pdb2gmx -f ../$(output_dir)/receptor_clean_previo.pdb -ignh -p ../$(output_dir)/receptor.top -i ../$(output_dir)/posre.itp -ff $(ff) -water $(water) -merge all
	gmx editconf -f conf.gro -o conf.pdb
	awk '{if(substr($$0,1,4)=="ATOM" || substr($$0,1,6)=="HETATM") $$0=substr($$0,1,21)"A"substr($$0,23); print}' conf.pdb > ../$(output_dir)/receptor.pdb
	python Scripts_Carlos/fix_pdb_chains.py ../$(output_dir)/receptor.pdb ../$(output_dir)/receptor_clean.pdb


	rm conf.gro conf.pdb

# make steps_previos_al_docking_ligando ligand_name=Abemaciclib root=/Users/danielcondetorres/Desktop/GUARDA_PROYECTO/MMDR/Code/Output_Carlos/experiment_33_docking_files/  resultados_vina_pdqt=3DKO_Abemaciclib_out.pdbqt.sdf
steps_previos_al_docking_ligando:
	cp $(root)/$(resultados_vina_pdqt) $(root)/$(ligand_name).sdf
	conda run -n acpype acpype -i  $(root)/$(ligand_name).sdf
	gmx editconf -f $(ligand_name).acpype/$(ligand_name)_GMX.gro -o ../$(output_dir)/md.pdb
#bash -c "cd $(basename $(notdir $(ligand_name))).acpype && source rungmx.sh"
#bash -c "cd $(basename $(notdir $(ligand_name))).acpype && gmx editconf -f md.gro -o md.pdb"

# make docking2D_resutls ligand_name=Abemaciclib output_dir=Prueba_Final root=/Users/danielcondetorres/Desktop/GUARDA_PROYECTO/MMDR/Code/Output_Carlos/experiment_33_docking_files/ receptor_pdb=3DKO.pdb
docking2D_resutls:
	mkdir ../$(output_dir)/Docking_picture/
	python Scripts_Carlos/docking2D.py --protein $(root)/$(receptor_pdb) --poses $(root)/$(ligand_name).sdf --output ../$(output_dir)/Docking_picture/interactions.html
#python Scripts_Carlos/docking2D.py --protein $(root)/receptor_clean.pdb --poses $(root)/$(ligand_name).sdf --output ../$(output_dir)/Docking_picture/interactions.html


# make combinar_topologias ligand_name=Abemaciclib output_dir=Prueba 
# make combinar_topologias ligand_name=Abemaciclib output_dir=Prueba 
combinar_topologias:
	cp $(ligand_name).acpype/$(ligand_name)_GMX.top ../$(output_dir)
	cp $(ligand_name).acpype/$(ligand_name)_GMX.itp ../$(output_dir)
	python Scripts_Carlos/combine_topols.py -l ../$(output_dir)/$(ligand_name)_GMX.top -r ../$(output_dir)/receptor.top  -o ../$(output_dir)/topol.top -ff $(ff) -itp $(ligand_name)

# make obtener_pdb_combinado  ligand_name=Abemaciclib output_dir=Prueba 
obtener_pdb_combinado:
	bash Scripts_Carlos/join_pdbqt2pdb.sh -r ../$(output_dir)/receptor_clean.pdb  -l $(ligand_name).acpype/$(ligand_name)_GMX.gro -o ../$(output_dir)/complejo.pdb

# make agua_e_iones output_dir=Prueba
agua_e_iones:
	gmx editconf -f ../$(output_dir)/complejo.pdb  -o ../$(output_dir)/box.pdb -d 1 # Necesito obtener de forma automatica las dimensiones de la caja
	#gmx solvate -cs spc216.gro -cp ../$(output_dir)/box.pdb -p ../$(output_dir)/topol.top -o ../$(output_dir)/caja_solvatada.pdb
	bash Scripts_Carlos/main_ligand.sh ../$(output_dir)/box.pdb ../$(output_dir)/	

analysis:
	bash Scripts_Carlos/analyisis_creation.sh ../$(output_dir)



# necesito todos los all.pdqt todas las conformacions para las gafas este flag dejar por ahora
all_pdqt:
	obabel -ipdbqt all.pdbqt -osdf -O all.sdf
	cp *pdbqt ../$(output_dir)
	cp *sdf ../$(output_dir)
	python Script_for_ligand_protein_formation/obtain_best_affinity.py -log log.txt
	cp log.txt ../$(output_dir)

#resultados_vina_pdqt sustituir por all.pdbqt 
# cambiar el python aa_interaction_flags.py por el nuevo de Alex, este debe estar viejo 
# make crear_archivos_VR ligand_name=Abemaciclib resultados_vina_pdqt=3DKO_Abemaciclib_out.pdbqt.sdf  output_dir=Prueba_final root=/Users/danielcondetorres/Desktop/GUARDA_PROYECTO/MMDR/Code/Output_Carlos/experiment_33_docking_files/ 
# Ojo cambiar por obabel -isdf all.pdqt.sdf -opdqt o lo que sea -O a.pdb segun los output que nos de esa gente en ligand_energy_attribution_bfcator.py
# cp $(root)/all.pdbqt ../$(output_dir)/all.pdbqt #Esta linea es la que vamos a necesitar, hacen falta los tests.

#cp  $(root)/$(resultados_vina_pdqt) ../$(output_dir)/all.pdbqt
crear_archivos_VR:
	cp $(root)/all.pdbqt ../$(output_dir)/
	./VR/bfactor.sh -p ../$(output_dir) -v ../$(output_dir)/VR -l ../$(output_dir)/md.pdb
	python VR/aa_interaction_dict.py --file ../$(output_dir) --pdb_path ../$(output_dir)/receptor_clean.pdb --output_path ../$(output_dir)/VR/aa_interactions_arx.json --sdf_path ../$(output_dir)/all.sdf
	python VR/aa_interaction_dict.py --file ../$(output_dir) --pdb_path ../$(output_dir)/$(receptor_pdb) --output_path ../$(output_dir)/VR/aa_interactions_arx_not_good.json --sdf_path ../$(output_dir)/all.sdf
	mkdir ../$(output_dir)/Simulation_files
	cp ../$(output_dir)/prod* ../$(output_dir)/Simulation_files/
# cp $(root)/all.pdbqt ../$(output_dir)/all.pdbqt #Esta linea es la que vamos a necesitar, hacen falta los tests.
clean:
	rm -r $(ligand_name).acpype *gro *mdp *edr *pdb *xtc *cpt *log
