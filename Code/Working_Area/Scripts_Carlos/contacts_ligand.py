import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import argparse
from itertools import combinations
from MDAnalysis.lib.distances import distance_array
import mdtraj as md

# ============================================================
# --------------------- HELPERS NUEVOS -----------------------
# ============================================================

def _choose_tick_step(n, target_labels=35):
    """
    Devuelve un paso entero para que haya ~target_labels etiquetas como máximo.
    Mantiene un escalado agradable (1,2,5,10,20,25,50,100,...).
    """
    if n <= 0:
        return 1
    if n <= target_labels:
        return 1
    nice_steps = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000]
    raw = max(1, int(np.ceil(n / target_labels)))
    step = next((s for s in nice_steps if s >= raw), nice_steps[-1])
    return step

def _set_time_ticks(ax, times, max_ticks=10):
    """
    Fija como mucho max_ticks etiquetas temporales, distribuidas uniformemente.
    """
    if len(times) == 0:
        return
    if len(times) <= max_ticks:
        ax.set_xticks(times)
    else:
        idx = np.linspace(0, len(times) - 1, num=max_ticks, dtype=int)
        ax.set_xticks([times[i] for i in idx])

def _place_compact_legend(ax, n_items):
    """
    Coloca la leyenda fuera del área del gráfico y en columnas si hay muchas curvas.
    """
    if n_items <= 8:
        ax.legend(fontsize=16, frameon=False)
    else:
        ncol = 2 if n_items <= 16 else 3
        ax.legend(
            fontsize=12,
            ncol=ncol,
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0,
            frameon=False
        )

# ============================================================
# ------------------- FUNCIONES EXISTENTES -------------------
# ============================================================

def calculate_contacts(pdb_file, xtc_file, contact_radius=3.5):
    # Cargar el sistema usando los archivos PDB y XTC
    u = mda.Universe(pdb_file, xtc_file, step=1)
    
    # Obtener la lista única de cadenas usando el atributo 'segid'
    chain_list = np.unique(u.atoms.segids)
    chain_list = [chain for chain in chain_list if chain and chain != "None"]
    
    # Preseleccionar los átomos correspondientes a cada cadena (se hace una sola vez)
    atom_groups = {chain: u.select_atoms(f"segid {chain}") for chain in chain_list}
    
    # Agregar el grupo SOL basado en resname SOL (si existe)
    sol_atoms = u.select_atoms("resname SOL")
    if len(sol_atoms) > 0:
        atom_groups["SOL"] = sol_atoms
    
    # Usar todas las claves del diccionario (cadenas + SOL)
    group_list = list(atom_groups.keys())
    
    # Inicializar un diccionario para almacenar los contactos entre cada par de grupos
    group_distances = {f"{g1}-{g2}": [] for g1, g2 in combinations(group_list, 2)}
    group_contacts = {f"{g1}-{g2}": [] for g1, g2 in combinations(group_list, 2)}
    times = []
    
    # Iterar sobre los fotogramas del archivo .xtc
    for ts in u.trajectory[::1]:
        print('Time: ', ts)
        times.append(ts.time / 1000.0)  # Convertir tiempo a ns (asumiendo que ts.time está en ps)
        for g1, g2 in combinations(group_list, 2):
            dists = distance_array(atom_groups[g1].positions, atom_groups[g2].positions, box=ts.dimensions)
            count = np.count_nonzero(dists < contact_radius)
            group_contacts[f"{g1}-{g2}"].append(count)
            avg_dist = np.mean(dists)  # Calcular la distancia promedio
            group_distances[f"{g1}-{g2}"].append(avg_dist)
    
    return group_contacts, times, group_distances

def plot_contacts(group_contacts, times, Output_dir, mapping):
    plt.figure(figsize=(12, 8))
    colors = plt.cm.get_cmap('tab20', len(group_contacts))
    
    for idx, (pair, contacts_list) in enumerate(group_contacts.items()):
        # Separar los identificadores de grupo
        g1, g2 = pair.split('-')
        # Aplicar el mapeo a cada grupo; si no se encuentra, se usa el valor original
        new_label = f"{mapping.get(g1, g1)}-{mapping.get(g2, g2)}"
        plt.plot(times, contacts_list, label=new_label, color=colors(idx))
    
    ax = plt.gca()
    _set_time_ticks(ax, times, max_ticks=10)                 # <<< NUEVO
    _place_compact_legend(ax, len(group_contacts))           # <<< NUEVO

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel("Time (ns)", fontsize=20)
    plt.ylabel("Number of Contacts", fontsize=20)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{Output_dir}/contactos_por_cadenas.png", bbox_inches="tight")
    print(f"Gráfico guardado como contactos_por_cadenas.png")
    plt.close()

def plot_distances(group_distances, times, Output_dir, mapping):
    plt.figure(figsize=(12, 8))
    colors = plt.cm.get_cmap('tab20', len(group_distances))
    
    for idx, (pair, distances_list) in enumerate(group_distances.items()):
        # Separar los identificadores de grupo
        g1, g2 = pair.split('-')
        # Aplicar el mapeo a cada grupo; si no se encuentra, se usa el valor original
        new_label = f"{mapping.get(g1, g1)}-{mapping.get(g2, g2)}"
        plt.plot(times, distances_list, label=new_label, color=colors(idx))
    
    ax = plt.gca()
    _set_time_ticks(ax, times, max_ticks=10)                 # <<< NUEVO
    _place_compact_legend(ax, len(group_distances))          # <<< NUEVO

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel("Time (ns)", fontsize=20)
    plt.ylabel("Average Distance (Å)", fontsize=20)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{Output_dir}/distances_por_cadenas.png", bbox_inches="tight")
    print(f"Gráfico guardado como distances_por_cadenas.png")
    plt.close()

def load_trajectory(pdb, xtc):
    """
    Load the system (structure and trajectory) into memory for faster subsequent analyses.

    Parameters:
    pdb (str): Path to the PDB file containing the protein structure.
    xtc (str): Path to the XTC file containing the trajectory.

    Returns:
    mda.Universe: Loaded MDAnalysis universe object containing both structure and trajectory.
    """
    u = mda.Universe(pdb, xtc)
    subunit1 = u.select_atoms("protein and (segid A or segid seg_0 or segid seg_0_Protein_chain_A or segid Protein_chain_A) and name CA")  # Cα atoms from subunit 1
    subunit3 =  u.select_atoms("segid C")
    if subunit3:  # False si vacío
        subunit2 = u.select_atoms("(segid B or segid seg_1 or segid seg_1_Protein_chain_B or segid Protein_chain_B or segid seg_1_Protein_chain_A2) and name CA")
    else:
        subunit2 = u.select_atoms("(segid B or segid seg_1 or segid seg_1_Protein_chain_B or segid Protein_chain_B or segid seg_1_Protein_chain_A2)")
    return u, subunit1, subunit2, subunit3

def contact_map(u, subunit_inter_1, subunit_inter_2, 
                res_init_1: int, res_init_2: int, 
                xlabel: str, ylabel: str, 
                Output_dir: str, contact_map_name: str, 
                cutoff: float = 5, jump_1: float = 3 , jump_2: float = 2,
                accordion: bool = True,         
                expand_factor: float = 2.0,     
                contact_threshold: float = 0.1,
                protected_zone: int = 5):       # zona de protección ±5 etiquetas
    """
    Calcula y visualiza el mapa de contactos entre dos subunidades.
    Si `accordion=True`, las filas con contacto > contact_threshold se dibujan más anchas,
    se muestran sus etiquetas, y se eliminan las etiquetas cercanas para evitar solapamiento.
    """

    for seg in u.segments:
        print(f"Segment (Unit): {seg.segid}, Name: {seg.residues.resnames}")

    residues1 = [f"{res.resname}{i+1+int(res_init_1)}" for i, res in enumerate(subunit_inter_1.residues)]
    try:
        residues2 = [f"{res.resname}{i+1+int(res_init_2)}" for i, res in enumerate(subunit_inter_2.residues)]
    except AttributeError:
        residues2 = []

    n_res_y = len(subunit_inter_2.residues)
    n_atoms_y = subunit_inter_2.atoms.n_atoms
    print(f'El NUMERO DE ATOMOS DEL RECEPTOR ES: {n_atoms_y}')

    # Grosor de línea dinámico
    linewidth_dict = {"lt_500": 0.5, "btw_500_1k": 0.2, "gt_1k": 0.08}
    if n_atoms_y > 1000:
        chosen_lw_key = "gt_1k"
        dynamic_jump_2 = 40
    elif n_atoms_y >= 500:
        chosen_lw_key = "btw_500_1k"
        dynamic_jump_2 = 10
    else:
        chosen_lw_key = "lt_500"
        dynamic_jump_2 = 5
    if n_res_y > 1000:
        dynamic_jump_2 = max(dynamic_jump_2, 20)
    chosen_lw = linewidth_dict[chosen_lw_key]
    print(f"[contact_map] Cadena eje Y -> residuos: {n_res_y}, linewidth elegido: {chosen_lw} ({chosen_lw_key})")

    # --- Calcular mapa de contactos ---
    contact_map_data = np.zeros((len(subunit_inter_2.atoms), len(subunit_inter_1.atoms)))
    normalizar = 0
    for ts in u.trajectory[::]:
        distances = distance_array(subunit_inter_2.positions, subunit_inter_1.positions)
        contacts = distances < cutoff
        contact_map_data += contacts.astype(int)
        normalizar += 1
    if normalizar > 0:
        contact_map_data /= normalizar

    # --- Determinar residuos con contacto ---
    has_contact = np.zeros(n_res_y, dtype=bool)
    for i, res in enumerate(subunit_inter_2.residues):
        atom_idx = np.where(subunit_inter_2.atoms.resindices == res.resindex)[0]
        if atom_idx.size > 0:
            if np.any(contact_map_data[atom_idx, :] > contact_threshold):
                has_contact[i] = True

    # --- Construir rejilla acordeón ---
    if accordion:
        weights = np.where(has_contact, expand_factor, 1.0)
        y_edges = np.concatenate(([0], np.cumsum(weights)))
    else:
        y_edges = np.arange(len(subunit_inter_2.atoms) + 1, dtype=float)
    x_edges = np.arange(len(subunit_inter_1.atoms) + 1, dtype=float)

    # --- Dibujar mapa ---
    fig, ax = plt.subplots(figsize=(18, 10))
    pcm = ax.pcolormesh(x_edges, y_edges, contact_map_data, cmap="plasma",
                        edgecolors="black", linewidth=chosen_lw, vmin=0, vmax=1)
    cbar = plt.colorbar(pcm, ax=ax, label="Contact Frequency")
    cbar.ax.tick_params(labelsize=14)
    ax.set_xlabel(f"{xlabel}", fontsize=18)
    ax.set_ylabel(f"{ylabel}", fontsize=18)

    # --- Etiquetas eje X ---
    ax.set_xticks(np.arange(0, len(residues1), jump_1))
    ax.set_xticklabels(residues1[::int(jump_1)], rotation=90, fontsize=11)

    # =======================
    # --- Etiquetas eje Y ---
    # =======================
    base_ticks = list(range(0, len(residues2), dynamic_jump_2))
    forced_ticks = list(np.where(has_contact)[0])
    all_ticks = sorted(set(base_ticks) | set(forced_ticks))

    if accordion:
        # Filtro de zona protegida alrededor de filas con contacto
        mask_ticks = np.ones(n_res_y, dtype=bool)
        for idx in forced_ticks:
            start = max(0, idx - protected_zone)
            end = min(n_res_y, idx + protected_zone + 1)
            mask_ticks[start:end] = False
        mask_ticks[forced_ticks] = True
        final_ticks = [i for i in all_ticks if mask_ticks[i]]

        # Posiciones de etiquetas (centros de filas)
        atom_resindices = subunit_inter_2.atoms.resindices
        first_atom_idx = [np.where(atom_resindices == res.resindex)[0][0] for res in subunit_inter_2.residues]
        y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
        y_ticks = [y_centers[i] for i in first_atom_idx if i < len(y_centers)]

        tick_positions = [y_ticks[i] for i in final_ticks if i < len(y_ticks)]
        tick_labels   = [residues2[i] for i in final_ticks if i < len(residues2)]
    else:
        # ✅ SIN ACORDEÓN: usa exactamente los mismos índices para posiciones y etiquetas
        final_ticks = base_ticks  # no mezclar con forced_ticks aquí para evitar desajustes
        tick_positions = final_ticks
        tick_labels   = [residues2[i] for i in final_ticks]

    ax.set_yticks(tick_positions)
    ax.set_yticklabels(tick_labels)

    # Formato etiquetas (solo aplica cuando hay acordeón y final_ticks marcado)
    if accordion:
        for i, lbl in zip(final_ticks, ax.get_yticklabels()):
            if has_contact[i]:
                lbl.set_fontweight("bold")
                lbl.set_color("black")
                lbl.set_fontsize(9)
            else:
                lbl.set_fontsize(6)
    else:
        for lbl in ax.get_yticklabels():
            lbl.set_fontsize(9)

    plt.tight_layout()
    plt.savefig(f"{Output_dir}/{contact_map_name}__.png", dpi=500)
    plt.close()

    return residues1, residues2


def dssp_chain_analysis(trajectory_file, topology_file, output_dir, name_file, chain_id, residues_names):
    """
    Analyze the secondary structure evolution of a specific protein chain over time using DSSP.
    """
    # Load trajectory and select only the specified chain
    traj = md.load(trajectory_file, top=topology_file)
    protein_atoms_indices = traj.topology.select(f'chainid {chain_id}')
    
    if len(protein_atoms_indices) == 0:
        print(f"Error: Chain {chain_id} not found in topology.")
        return

    traj_protein_slice = traj.atom_slice(protein_atoms_indices) 

    # Convert time from picoseconds to nanoseconds
    time_ns = traj.time / 1000  # Convert ps to ns

    # Compute DSSP for the selected chain
    dssp = md.compute_dssp(traj_protein_slice)

    # DSSP mapping
    dssp_mapping = {
        'H': 0, 'G': 0, 'I': 0,  # Helix
        'E': 1, 'B': 1,          # Beta sheet
        'T': 2, 'S': 2, 'C': 2   # Coil
    }

    # Convert DSSP characters to numeric matrix
    dssp_numeric = np.vectorize(lambda x: dssp_mapping.get(x, 2))(dssp)

    # Plot heatmap
    plt.figure(figsize=(14, 6))
    plt.imshow(dssp_numeric.T, aspect='auto', cmap="coolwarm", interpolation='nearest')

    # Labels
    plt.xlabel("Time (ns)", fontsize=18)
    plt.ylabel("Residues", fontsize=18)

    # Eje X (tiempo) con número limitado de marcas
    frame_indices = np.linspace(0, len(time_ns) - 1, num=10, dtype=int)
    plt.xticks(frame_indices, np.round(time_ns[frame_indices], 1), fontsize=14)

    # Eje Y (residuos) con salto dinámico uniforme (mismo criterio que en el resto)
    n_res = len(residues_names)
    step = _choose_tick_step(n_res, target_labels=35)
    yticks_idx = np.arange(0, n_res, step) if n_res > 0 else np.array([0])
    if n_res > 0:
        plt.yticks(yticks_idx, np.array(residues_names)[yticks_idx], rotation=0, fontsize=12)

    # Colorbar y leyenda
    cbar = plt.colorbar(label="Secondary Structure")
    cbar.set_ticks([0, 1, 2])
    cbar.set_ticklabels(["Helix", "Beta sheet", "Coil"])
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.yaxis.label.set_size(16)

    # Guardar
    output_path = f"{output_dir}/{name_file}_chain_{chain_id}_dssp.png"
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()  # Close figure to avoid overlap

# ============================================================
# --------------------------- MAIN ---------------------------
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Calcular el número de contactos entre cadenas (segids) y con SOL a lo largo del tiempo en una simulación."
    )
    parser.add_argument('-pdb', type=str, required=True, help="Archivo PDB con la topología")
    parser.add_argument('-xtc', type=str, required=True, help="Archivo XTC con la trayectoria")
    parser.add_argument('--contact_radius', type=float, default=8, 
                        help="Radio de contacto en Ångströms (por defecto 3.5 Å)")
    parser.add_argument('-Output_dir', type=str, default='.', 
                        help="Nombre del archivo de salida para el gráfico")
    parser.add_argument("-r1_number", default=1, help="Int number to start to count residue 1")
    parser.add_argument("-r2_number", default=1, help="Int number to start to count residue 2")
    parser.add_argument("-r3_number", default=0, help="Int number to start to count residue 3")
    parser.add_argument("-label_1", default='Protein', help="Name residue 1")
    parser.add_argument("-label_2", default='Ligand', help="Name residue 2")
    parser.add_argument("-label_3", default='Other thing', help="Name residue 3")
    args = parser.parse_args()

    # Definir el mapeo de segids a los nombres deseados.
    mapping = {"A": args.label_1, "B": args.label_2}

    group_contacts, times, group_distances = calculate_contacts(args.pdb, args.xtc, args.contact_radius)
    plot_contacts(group_contacts, times, args.Output_dir, mapping)
    plot_distances(group_distances, times, args.Output_dir, mapping)

    u, s1, s2, s3 = load_trajectory(args.pdb, args.xtc)

    residues_1_names, residues_2_names = contact_map(
        u, s2, s1,
        args.r2_number, args.r1_number,
        args.label_2, args.label_1,
        args.Output_dir, '1',
        args.contact_radius,
        jump_1=10,
        jump_2=10,
        accordion=True,           # activa el acordeón
        expand_factor=30,         # factor de expansión
        contact_threshold=0.01    # umbral mínimo para “hay contacto”
    )

    residues_1_names, residues_2_names = contact_map(
        u, s2, s1,
        args.r2_number, args.r1_number,
        args.label_2, args.label_1,
        args.Output_dir, 'not_accordion',
        args.contact_radius,
        jump_1=10,
        jump_2=10,
        accordion=False,           # activa el acordeón
        expand_factor=20,         # factor de expansión
        contact_threshold=0.01    # umbral mínimo para “hay contacto”
    )

    dssp_chain_analysis(args.xtc, args.pdb, args.Output_dir, args.label_1, 0, residues_2_names)

if __name__ == "__main__":
    main()
