import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import argparse
from itertools import combinations
from MDAnalysis.lib.distances import distance_array
import mdtraj as md
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
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel("Time (ns)", fontsize=20)
    plt.ylabel("Number of Contacts", fontsize=20)
    plt.legend(fontsize=20)
    plt.grid(True)
    plt.savefig(f"{Output_dir}/contactos_por_cadenas.png")
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
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel("Time (ns)", fontsize=20)
    plt.ylabel("Average Distance (Å)", fontsize=20)
    plt.legend(fontsize=20)
    plt.grid(True)
    plt.savefig(f"{Output_dir}/distances_por_cadenas.png")
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

    subunit3 =  u.select_atoms("segid C")#.center_of_mass()
    if subunit3:  # False if empty
        subunit2 = u.select_atoms("(segid B or segid seg_1 or segid seg_1_Protein_chain_B or segid Protein_chain_B or segid seg_1_Protein_chain_A2) and name CA")
    # If not subunit3 it means that subunit 2 is the ligand, so we dont need CA
    else:
        subunit2 = u.select_atoms("(segid B or segid seg_1 or segid seg_1_Protein_chain_B or segid Protein_chain_B or segid seg_1_Protein_chain_A2)")
    return u, subunit1, subunit2, subunit3


def contact_map(u, subunit_inter_1, subunit_inter_2, res_init_1: int, res_init_2: int, xlabel: str, ylabel: str, Output_dir: str, contact_map_name: str, cutoff: float = 5, jump_1: float = 3 , jump_2: float = 2):
    """
    This function calculates and visualizes the contact map between two protein subunits 
    based on their Cα atoms from a trajectory.

    Parameters:
    u (mda.Universe): Loaded MDAnalysis universe object containing both structure and trajectory.
    cutoff (float): Distance threshold (in Angstroms) for considering a contact. Default is 12 Å.
    frames_to_analyze (int): Number of frames to analyze from the trajectory. Default is 1 frame.

    Returns:
    None: Displays a contact map as a heatmap.
    """

    # Print information about segments (molecules/chains)
    for seg in u.segments:
        print(f"Segment (Unit): {seg.segid}, Name: {seg.residues.resnames}")

    # Select subunits (adjust the segment IDs based on your system)
    # Get the residue names and renumber them starting from 1
    residues1 = [f"{res.resname}{i+1+int(res_init_1)}" for i, res in enumerate(subunit_inter_1.residues)]
    try:
        residues2 = [f"{res.resname}{i+1+int(res_init_2)}" for i, res in enumerate(subunit_inter_2.residues)]
        print('Res', residues2)
    except AttributeError:
        pass
    # Initialize the contact map (matrix)
    contact_map = np.zeros((len(subunit_inter_2.atoms), len(subunit_inter_1.atoms)))
    normalizar = 0
    # Calculate contacts for each frame
    for ts in u.trajectory[::]:
        # Calculate the distance array between residues of subunit 1 and subunit 2
        distances = distance_array(subunit_inter_2.positions, subunit_inter_1.positions)
        #print('Distances: ',subunit_inter_1.positions)
        # Define contacts as those with distance less than the cutoff
        contacts = distances < cutoff
        #print(contacts)
        # Update the contact map (sum contacts in each frame)
        contact_map += contacts.astype(int)
        normalizar += 1
    # Normalize the contact map
    contact_map /= normalizar

    # Create the contact map figure
    plt.figure(figsize=(18, 10))
    plt.pcolormesh(contact_map, cmap="plasma", edgecolors="black", linewidth=0.5,vmin=0, vmax=1)

    # Add color bar for contact frequency
    cbar = plt.colorbar(label="Contact Frequency")
    cbar.ax.tick_params(labelsize=14)

    # Force all residues to be labeled on the axes
    plt.xticks(ticks=np.arange(len(residues1)), labels=residues1, rotation=90, fontsize=8)
    plt.yticks(ticks=np.arange(len(residues2)), labels=residues2, fontsize=8)
    plt.gca().set_xticks(np.arange(len(residues1)))  # Ensure every tick has a label
    plt.gca().set_yticks(np.arange(len(residues2)))

    # Labels and title
    plt.xlabel(f"{xlabel}", fontsize=18)
    plt.ylabel(f"{ylabel}", fontsize=18)
        # Adjust xticks and yticks spacing to make the plot more readable
    plt.gca().set_xticks(np.arange(len(residues1)))  # Ensure every tick has a label
    plt.gca().set_yticks(np.arange(len(residues2)))
    plt.xticks(ticks=np.arange(0, len(residues1), jump_1), labels=residues1[::jump_1], rotation=90, fontsize=11)  # Show every 5th label
    plt.yticks(ticks=np.arange(0, len(residues2), jump_2), labels=residues2[::jump_2], fontsize=11)  # Show every 5th label

    plt.tight_layout()
    # Show the plot
    plt.savefig(f"{Output_dir}/{contact_map_name}__.png", dpi = 500)
    plt.close()
    return residues1, residues2

def dssp_chain_analysis(trajectory_file, topology_file, output_dir, name_file, chain_id,residues_names):
    """
    Analyze the secondary structure evolution of a specific protein chain over time using DSSP.

    Parameters:
    trajectory_file (str): Path to the trajectory file (.xtc, .dcd, etc.).
    topology_file (str): Path to the topology file (.pdb, .gro, etc.).
    output_dir (str): Directory to save the plot.
    name_file (str): Name for the output file.
    chain_id (int): Chain ID to analyze.

    Returns:
    None: Saves a heatmap of the secondary structure evolution for the specified chain.
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

    # Adjust x-axis ticks to show time in ns
    frame_indices = np.linspace(0, len(time_ns) - 1, num=10, dtype=int)
    plt.xticks(frame_indices, np.round(time_ns[frame_indices], 1),fontsize = 14)

    # Colorbar and legend
    cbar = plt.colorbar(label="Secondary Structure")
    cbar.set_ticks([0, 1, 2])
    cbar.set_ticklabels(["Helix", "Beta sheet", "Coil"])
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.yaxis.label.set_size(16)  # Ajusta el tamaño del label
    # Save plot
    output_path = f"{output_dir}/{name_file}_chain_{chain_id}_dssp.png"
    plt.yticks(ticks=np.arange(0, len(residues_names), 20), labels=residues_names[::20], rotation=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()  # Close figure to avoid overlap



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
    # Por ejemplo, si se detectan segids A, B y C, se asignan los nombres deseados.
    #mapping = {"A": args.label_1, "B": args.label_2, "C": args.label_3}

    mapping = {"A": args.label_1, "B": args.label_2}
    group_contacts, times, group_distances = calculate_contacts(args.pdb, args.xtc, args.contact_radius)
    plot_contacts(group_contacts, times, args.Output_dir, mapping)
    plot_distances(group_distances, times,  args.Output_dir, mapping)
    u, s1, s2, s3 = load_trajectory(args.pdb, args.xtc)
    residues_1_names, residues_2_names = contact_map(u, s2, s1, args.r2_number, args.r1_number, args.label_2, args.label_1, args.Output_dir, '1', args.contact_radius,10,10)
    #try:
    #    contact_map(u, s2, s3, args.r2_number, args.r3_number, args.label_2, args.label_3, args.Output_dir,'2', args.contact_radius,1 ,1)
    #    contact_map(u, s1, s3, args.r1_number, args.r3_number, args.label_1, args.label_3, args.Output_dir,'3', args.contact_radius,4 ,1)
    #except:
    #    print('Only one protein')
    dssp_chain_analysis(args.xtc, args.pdb, args.Output_dir, args.label_1, 0, residues_2_names)
if __name__ == "__main__":
    main()
