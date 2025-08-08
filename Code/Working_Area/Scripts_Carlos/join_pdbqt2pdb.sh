#!/bin/bash

# Valores por defecto
output="complex.pdb"

# Leer flags
while getopts r:l:o: flag; do
  case "${flag}" in
    r) receptor_pdbqt=${OPTARG};;
    l) ligand_pdbqt=${OPTARG};;
    o) output=${OPTARG};;
  esac
done

# Verificar inputs
if [[ -z "$receptor_pdbqt" || -z "$ligand_pdbqt" ]]; then
  echo "❌ Uso: $0 -r receptor.pdbqt -l ligando.pdbqt -o salida.pdb"
  exit 1
fi

# Convertir a PDB
obabel "$receptor_pdbqt" -O temp_receptor.pdb || { echo "❌ Error convirtiendo receptor"; exit 1; }
obabel "$ligand_pdbqt" -O temp_ligand.pdb || { echo "❌ Error convirtiendo ligando"; exit 1; }

# Limpiar líneas END
grep -v -e '^END' -e '^TER' temp_receptor.pdb > receptor_clean.pdb
grep -v -e '^END' -e '^TER' temp_ligand.pdb > ligand_clean.pdb

# Combinar
cat receptor_clean.pdb ligand_clean.pdb > "$output"
echo "END" >> "$output"

# Limpiar
rm temp_receptor.pdb temp_ligand.pdb receptor_clean.pdb ligand_clean.pdb

echo "✅ Archivo combinado generado: $output"

