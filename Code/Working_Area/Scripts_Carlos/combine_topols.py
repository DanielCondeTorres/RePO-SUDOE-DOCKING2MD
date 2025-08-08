import argparse

def extract_molecules(lines):
    """Extrae moléculas de la sección [ molecules ]"""
    in_molecules = False
    molecules = []
    for line in lines:
        if line.strip().startswith('[ molecules ]'):
            in_molecules = True
            continue
        if in_molecules:
            if line.strip() == '' or line.strip().startswith('['):
                break
            if not line.strip().startswith(';'):
                molecules.append(line.strip())
    return molecules

def merge_topologies(receptor_top, ligand_top, output_top):
    with open(receptor_top, 'r') as f:
        receptor_lines = f.readlines()

    with open(ligand_top, 'r') as f:
        ligand_lines = f.readlines()

    output_lines = []

    abem_itp_line = '#include "Abemaciclib_GMX.itp"\n'
    posre_line = '#include "posre.itp"\n'

    # Etapas de control
    abem_itp_inserted = False
    posre_inserted = False
    inside_posres_block = False
    skip_rest = False

    for i, line in enumerate(receptor_lines):
        stripped = line.strip()

        if stripped.startswith('[ system ]') or stripped.startswith('[ molecules ]'):
            skip_rest = True
            break  # Salta todo lo que sigue después

        # 1. Insert Abemaciclib_GMX.itp justo después de forcefield.itp
        if stripped == '#include "charmm27.ff/forcefield.itp"' and not abem_itp_inserted:
            output_lines.append(line)
            output_lines.append(abem_itp_line)
            abem_itp_inserted = True
            continue

        # 2. Insert posre.itp dentro del bloque POSRES
        if stripped == '#ifdef POSRES':
            output_lines.append(line)
            inside_posres_block = True
            continue
        elif stripped == '#endif' and inside_posres_block:
            if not posre_inserted:
                output_lines.append(posre_line)
                posre_inserted = True
            output_lines.append(line)
            inside_posres_block = False
            continue

        output_lines.append(line)

    output_lines.append("\n")

    # 3. Añadir [ system ]
    output_lines.append('[ system ]\n')
    output_lines.append(' Protein-Ligand Complex\n\n')

    # 4. Añadir [ molecules ]
    output_lines.append('[ molecules ]\n')

    receptor_mols = extract_molecules(receptor_lines)
    ligand_mols = extract_molecules(ligand_lines)

    seen = set()
    for mol in receptor_mols + ligand_mols:
        name = mol.split()[0]
        if name not in seen:
            output_lines.append(f"{mol}\n")
            seen.add(name)

    # 5. Escribir archivo final
    with open(output_top, 'w') as f:
        f.writelines(output_lines)

    print(f"✅ Archivo combinado creado: {output_top}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combina topologías GROMACS de receptor y ligando.")
    parser.add_argument("-r", "--receptor", required=True, help="Archivo .top del receptor")
    parser.add_argument("-l", "--ligand", required=True, help="Archivo .top del ligando")
    parser.add_argument("-o", "--output", required=True, help="Archivo .top de salida combinado")

    args = parser.parse_args()
    merge_topologies(args.receptor, args.ligand, args.output)


'''
def extract_molecules(lines):
    """Extrae moléculas de la sección [ molecules ]"""
    in_molecules = False
    molecules = []
    for line in lines:
        if line.strip().startswith('[ molecules ]'):
            in_molecules = True
            continue
        if in_molecules:
            if line.strip() == '' or line.strip().startswith('['):
                break
            if not line.strip().startswith(';'):
                molecules.append(line.strip())
    return molecules

def merge_topologies(receptor_top, ligand_top, output_top):
    with open(receptor_top, 'r') as f:
        receptor_lines = f.readlines()

    with open(ligand_top, 'r') as f:
        ligand_lines = f.readlines()

    output_lines = []

    # 1. [ defaults ] from ligand
    copying = False
    for line in ligand_lines:
        if line.strip().startswith('[ defaults ]'):
            copying = False
        if copying:
            output_lines.append(line)
            if line.strip() == '':
                break

    output_lines.append("\n")

    # 2. All forcefield and include lines from receptor
    includes_set = set()
    for line in receptor_lines:
        if line.strip().startswith('#include'):
            if line not in includes_set:
                output_lines.append(line)
                includes_set.add(line)

    output_lines.append("\n")

    # 3. Content before [ system ] in receptor (e.g., moleculetype, atoms, etc.)
    for line in receptor_lines:
        if line.strip().startswith('[ system ]'):
            break
        if not line.strip().startswith('#include'):
            output_lines.append(line)

    output_lines.append("\n")

    # 4. Add ligand .itp include
    for line in ligand_lines:
        if 'Abemaciclib_GMX.itp' in line and line.strip().startswith('#include'):
            output_lines.append(line)
            break

    # 5. Add posre include if present
    for line in ligand_lines:
        if 'posre_Abemaciclib.itp' in line and line.strip().startswith('#include'):
            output_lines.append(line)
            break

    output_lines.append("\n")

    # 6. Write [ system ]
    output_lines.append('[ system ]\n')
    output_lines.append(' Protein-Ligand Complex\n\n')

    # 7. Write [ molecules ]
    output_lines.append('[ molecules ]\n')

    # Add unique molecule entries: receptor first, ligand second
    receptor_mols = extract_molecules(receptor_lines)
    ligand_mols = extract_molecules(ligand_lines)

    seen = set()
    for mol in receptor_mols + ligand_mols:
        name = mol.split()[0]
        if name not in seen:
            output_lines.append(f"{mol}\n")
            seen.add(name)

    # 8. Write output file
    with open(output_top, 'w') as f:
        f.writelines(output_lines)

    print(f"✅ Archivo combinado creado: {output_top}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combina topologías GROMACS de receptor y ligando.")
    parser.add_argument("-r", "--receptor", required=True, help="Archivo .top del receptor")
    parser.add_argument("-l", "--ligand", required=True, help="Archivo .top del ligando")
    parser.add_argument("-o", "--output", required=True, help="Archivo .top de salida combinado")

    args = parser.parse_args()
    merge_topologies(args.receptor, args.ligand, args.output)

'''