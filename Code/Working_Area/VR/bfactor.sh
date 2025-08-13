#!/bin/bash

# ================================
# Parse command-line arguments
# ================================
while getopts p:v:l: flag
do
    case "${flag}" in
        p) path=${OPTARG};;
        v) vr_folder=${OPTARG};;
        l) ligand=${OPTARG};;
        *) echo "Usage: $0 -p <path> -v <vr_folder> -l <ligand>"
           exit 1;;
    esac
done

if [ -z "$path" ] || [ -z "$vr_folder" ] || [ -z "$ligand" ]; then
    echo "Missing required arguments."
    echo "Usage: $0 -p <path> -v <vr_folder> -l <ligand>"
    exit 1
fi

# ================================
# Activate Conda environment
# ================================

# Create output folder
mkdir -p "$vr_folder"

# Convert and concatenate trajectories
gmx convert-trj -f "$path/min_fep1.trr" -o "$path/min_fep1_c.xtc"

gmx trjcat -f "$path/min_fep1_c.xtc" "$path/eq_nvt_fep.xtc" "$path/eq_npt_fep_1.xtc" "$path/prod.xtc" -o "$path/traj.xtc"
gmx make_ndx -f "$path/prod.tpr" -o "$path/index.ndx" << EOF
 "Protein" | "Other"
q
EOF
gmx trjconv -s "$path/prod.tpr" -f "$path/traj.xtc" -pbc mol -center -n "$path/index.ndx" -o "$path/center_complete.pdb" -e 0 << EOF
Protein 
Protein_Other
q
EOF

gmx trjconv -s "$path/prod.tpr" -f "$path/traj.xtc" -pbc mol -center -n "$path/index.ndx" -o "$path/center_complete.xtc" -b 0 << EOF
Protein 
Protein_Other
q
EOF

# Execute the Python script
python VR/ligand_energy_attribution_bfactor.py \
    --pdbqt_ligs  "$path/all.pdbqt" \
    --output_pdb "$path/all" \
    --receptor "$path/center_complete.pdb" \
    --dir_models "$path/models_center_pdb" \
    --dir_output "$path/bfactor_receptor_center" \
    --cutoff 5.0

# Copy and rename output files
cp "$path/all.pdb" "$vr_folder"
cp "$path/bfactor_receptor_center/receptor_final_bfactor_residue.pdb" "$vr_folder"
cp "$path/center_complete.xtc" "$vr_folder/center.xtc"
mv "$vr_folder/receptor_final_bfactor_residue.pdb" "$vr_folder/center.pdb"

# Copy ligand pdb file
cp $ligand "$vr_folder/ligando.pdb"
