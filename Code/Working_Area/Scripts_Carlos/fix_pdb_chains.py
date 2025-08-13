# fix_pdb_resids.py
import sys

if len(sys.argv) != 3:
    print(f"Uso: python {sys.argv[0]} input.pdb output.pdb")
    sys.exit(1)

input_pdb = sys.argv[1]
output_pdb = sys.argv[2]

new_lines = []
resid_offset = 0
prev_resid = None
max_resid_seen = 0

with open(input_pdb, 'r') as f:
    for line in f:
        if line.startswith(("ATOM", "HETATM")):
            resid_str = line[22:26]  # columnas resid
            try:
                resid = int(resid_str)
            except ValueError:
                new_lines.append(line)
                continue

            # Detectar reinicio
            if prev_resid is not None and resid < prev_resid:
                resid_offset = max_resid_seen
            new_resid = resid + resid_offset

            # Guardar máximo visto
            if new_resid > max_resid_seen:
                max_resid_seen = new_resid

            # Reemplazar en la línea (resid ocupa columnas 23–26, ancho 4)
            new_line = line[:22] + f"{new_resid:4d}" + line[26:]
            new_lines.append(new_line)
            prev_resid = resid
        else:
            new_lines.append(line)

with open(output_pdb, 'w') as f:
    f.writelines(new_lines)

print(f"PDB corregido guardado en {output_pdb}")

