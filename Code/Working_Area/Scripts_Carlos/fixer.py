#!/usr/bin/env python3

import argparse
import os
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def clean_pdb_for_gromacs(input_file, output_file, keep_waters=False, keep_ions=False, add_missing_residues=False):
    """
    Limpia y repara archivo PDB para uso con GROMACS
    """
    print(f"🔧 Procesando {input_file}...")
    
    # Inicializar PDBFixer
    try:
        fixer = PDBFixer(filename=input_file)
    except Exception as e:
        print(f"❌ Error al leer el archivo PDB: {e}")
        sys.exit(1)
    
    # Mostrar información inicial
    print(f"📊 Cadenas encontradas: {[chain.id for chain in fixer.topology.chains()]}")
    
    # Reparar estructura
    print("🔍 Buscando residuos faltantes...")
    fixer.findMissingResidues()
    
    if not add_missing_residues:
        print("⚠️  Evitando añadir residuos faltantes (para evitar cadenas extrañas)")
        fixer.missingResidues = {}  # Esto evita las cadenas extrañas
    else:
        print("➕ Añadiendo residuos faltantes (esto puede crear cadenas largas)")
    
    print("🔍 Buscando átomos faltantes...")
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    
    print("💧 Añadiendo hidrógenos (pH=7.0)...")
    fixer.addMissingHydrogens(pH=7.0)
    
    # Guardar archivo temporal
    temp_pdb = output_file + ".tmp"
    print("💾 Guardando archivo temporal...")
    
    try:
        with open(temp_pdb, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
    except Exception as e:
        print(f"❌ Error al escribir archivo temporal: {e}")
        sys.exit(1)
    
    # Limpiar archivo: eliminar HETATM no deseados
    print("🧹 Limpiando HETATM no deseados...")
    
    # Iones comunes que podrías querer mantener
    ions_to_keep = {'NA', 'CL', 'MG', 'CA', 'K', 'ZN', 'FE'} if keep_ions else set()
    
    with open(temp_pdb, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith("HETATM"):
                # Extraer nombre del residuo (posiciones 17-20)
                res_name = line[17:20].strip()
                
                # Mantener aguas si se especifica
                if keep_waters and res_name in ['HOH', 'WAT', 'TIP3']:
                    f_out.write(line)
                # Mantener iones si se especifica
                elif keep_ions and res_name in ions_to_keep:
                    f_out.write(line)
                # Eliminar otros HETATM
                else:
                    continue
            else:
                # Mantener todas las líneas que no sean HETATM
                f_out.write(line)
    
    # Limpiar archivo temporal
    os.remove(temp_pdb)
    
    print(f"✅ Archivo limpio guardado como: {output_file}")
    
    # Mostrar estadísticas del archivo final
    show_pdb_stats(output_file)

def show_pdb_stats(pdb_file):
    """Muestra estadísticas del archivo PDB procesado"""
    atom_count = 0
    chain_ids = set()
    res_types = set()
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom_count += 1
                chain_ids.add(line[21])
                res_types.add(line[17:20].strip())
    
    print(f"\n📈 Estadísticas del archivo final:")
    print(f"   • Átomos totales: {atom_count}")
    print(f"   • Cadenas: {sorted(chain_ids)}")
    print(f"   • Tipos de residuos: {len(res_types)}")

def main():
    parser = argparse.ArgumentParser(
        description='Repara y limpia archivos PDB para uso con GROMACS',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
  python pdb_cleaner.py -i protein.pdb -o protein_clean.pdb
  python pdb_cleaner.py -i protein.pdb -o protein_clean.pdb --keep-waters --keep-ions
        """
    )
    
    parser.add_argument('-i', '--input', required=True, 
                       help='Archivo PDB de entrada')
    parser.add_argument('-o', '--output', required=True, 
                       help='Archivo PDB de salida limpio')
    parser.add_argument('--add-missing-residues', action='store_true',
                       help='Añadir residuos faltantes (puede crear cadenas largas)')
    parser.add_argument('--keep-waters', action='store_true',
                       help='Mantener moléculas de agua (HOH, WAT)')
    parser.add_argument('--keep-ions', action='store_true',
                       help='Mantener iones comunes (Na+, Cl-, Mg2+, etc.)')
    parser.add_argument('--no-missing-residues', action='store_true',
                       help='No añadir residuos faltantes al inicio/final (DEPRECADO - usar por defecto)')
    
    args = parser.parse_args()
    
    # Verificar que el archivo de entrada existe
    if not os.path.exists(args.input):
        print(f"❌ Error: No se encuentra el archivo {args.input}")
        sys.exit(1)
    
    # Procesar archivo
    clean_pdb_for_gromacs(
        args.input, 
        args.output, 
        keep_waters=args.keep_waters,
        keep_ions=args.keep_ions,
        add_missing_residues=args.add_missing_residues
    )
    
    print(f"\n🎉 ¡Proceso completado! Ahora puedes usar:")
    print(f"   gmx pdb2gmx -f {args.output} -o processed.gro -p topol.top")

if __name__ == "__main__":
    main()



'''
import argparse
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def main():
    parser = argparse.ArgumentParser(description='Repara PDB y elimina HETATM para GROMACS.')
    parser.add_argument('-i', '--input', required=True, help='Archivo PDB de entrada')
    parser.add_argument('-o', '--output', required=True, help='Archivo PDB de salida limpio')
    args = parser.parse_args()

    fixer = PDBFixer(filename=args.input)

    # Repara la estructura
    fixer.findMissingResidues()
    #fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    # Guarda archivo temporal
    temp_pdb = args.output + ".tmp"
    with open(temp_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    # Ahora abre el temp, elimina todas las líneas HETATM y guarda final
    with open(temp_pdb, 'r') as f_in, open(args.output, 'w') as f_out:
        for line in f_in:
            if not line.startswith("HETATM"):
                f_out.write(line)

    print(f"✅ Archivo reparado y limpio guardado como: {args.output}")

if __name__ == "__main__":
    main()

'''