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

