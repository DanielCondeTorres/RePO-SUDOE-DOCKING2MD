#!/bin/bash

OUTPUT_DIR="$1"  # Archivo PDB del complejo
module load cesga/2020 gcc/system openmpi/4.0.5_ft3_cuda gromacs/2021.4-plumed-2.8.0

gmx make_ndx -f "$OUTPUT_DIR"/prod.tpr -o "$OUTPUT_DIR"/index.ndx << EOF
"Other"|"Protein"
 q
EOF
gmx trjconv -s "$OUTPUT_DIR"/prod.tpr -f "$OUTPUT_DIR"/prod.xtc -pbc mol -center -n "$OUTPUT_DIR"/index.ndx -o "$OUTPUT_DIR"/center.pdb -e 0 << EOF
Protein
Other_Protein
EOF
gmx trjconv -s "$OUTPUT_DIR"/prod.tpr -f "$OUTPUT_DIR"/prod.xtc -pbc mol -center -n "$OUTPUT_DIR"/index.ndx -o "$OUTPUT_DIR"/center.xtc -skip 10 << EOF
Protein
Other_Protein
EOF
mkdir "$OUTPUT_DIR"/Images/

module load cesga/system miniconda3/22.11.1-1
conda run -n acpype python ../Working_Area/Scripts_Carlos/analysis_complex.py -pdb "$OUTPUT_DIR"/center.pdb -xtc "$OUTPUT_DIR"/center.xtc -output "$OUTPUT_DIR"/Images/
conda run -n acpype python ../Working_Area/Scripts_Carlos/contacts_ligand.py -pdb "$OUTPUT_DIR"/center.pdb -xtc "$OUTPUT_DIR"/center.xtc -Output_dir "$OUTPUT_DIR"/Images/


