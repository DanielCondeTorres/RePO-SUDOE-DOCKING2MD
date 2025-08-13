import numpy as np
import matplotlib.pyplot as plt
import argparse
from MDAnalysis.analysis import rms
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import mdtraj as md
import matplotlib.patches as mpatches
import matplotlib as mpl
from MDAnalysis.analysis import align
mpl.use('Agg')

def moving_average(x: list, n: int=30):
    """
    This function computes the moving average of a time series (array `x`) over a specified window size `n`.

    Parameters:
    x (array-like): Input time series data. It is a list or NumPy array of numerical values.
    n (int): The window size for the moving average. This is the number of consecutive data points used to calculate the average.

    Returns:
    list: A list of moving average values, where each value represents the average of the `n` preceding data points in the time series.

    Notes:
    - The first `n-1` values will be `NaN` or undefined, since there are not enough points at the beginning to calculate a full window average.
    - This method uses a cumulative sum approach to efficiently compute the moving average.
    """

    # Initialize the result list
    result = []
    # Calculate the moving average for the first n-1 points (undefined)
    for i in range(n):
        if i == 0:
            pass  # The first element has no moving average, so we skip it
        else:
            # Calculate the average of the first i points
            result.append(float(np.sum(x[0:i]) / len(x[0:i])))
    # Calculate cumulative sum of the input array
    cumsum = np.cumsum(np.insert(x, 0, 0))
    # Apply the formula to compute the moving average
    res = (cumsum[n:] - cumsum[:-n]) / float(n)
    # Append the result of the moving average to the result list
    for elem in res:
        result.append(elem)
    return result




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
    return u

def contact_map(u, res_init_1: int, res_init_2: int, xlabel: str, ylabel: str,Output_dir: str, cutoff: float = 5, frames_to_analyze: int = 1):
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

    # Access the trajectory frames
    for ts in u.trajectory:  # Iterate over the first 10 frames
        print(f"Frame {ts.frame}: Time {ts.time} ps")
    print('RES: ',res_init_1)
    # Select subunits (adjust the segment IDs based on your system)
    subunit1 = u.select_atoms("protein and (segid A or segid seg_0 or segid seg_0_Protein_chain_A or segid Protein_chain_A) and name CA")  # Cα atoms from subunit 1
    print('subunit1:', subunit1)
    subunit2 = u.select_atoms("protein and (segid B or segid seg_1 or segid seg_1_Protein_chain_B or segid Protein_chain_B or segid seg_1_Protein_chain_A2) and name CA")  # Cα atoms from subunit 2

    # Get the residue names and renumber them starting from 1
    residues1 = [f"{res.resname}{i+1+int(res_init_1)}" for i, res in enumerate(subunit1.residues)]
    residues2 = [f"{res.resname}{i+1+int(res_init_2)}" for i, res in enumerate(subunit2.residues)]

    # Initialize the contact map (matrix)
    contact_map = np.zeros((len(subunit2.residues), len(subunit1.residues)))
    normalizar = 0
    # Calculate contacts for each frame
    for ts in u.trajectory[::]:
        # Calculate the distance array between residues of subunit 1 and subunit 2
        distances = distance_array(subunit2.positions, subunit1.positions)
        
        # Define contacts as those with distance less than the cutoff
        contacts = distances < cutoff
        
        # Update the contact map (sum contacts in each frame)
        contact_map += contacts.astype(int)
        normalizar += 1
    # Normalize the contact map
    contact_map /= normalizar

    # Create the contact map figure
    plt.figure(figsize=(18, 10))
    plt.pcolormesh(contact_map, cmap="plasma", edgecolors="black", linewidth=0.5)

    # Add color bar for contact frequency
    plt.colorbar(label="Contact Frequency")

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
    plt.xticks(ticks=np.arange(0, len(residues1), 4), labels=residues1[::4], rotation=90, fontsize=11)  # Show every 5th label
    plt.yticks(ticks=np.arange(0, len(residues2), 1), labels=residues2[::1], fontsize=11)  # Show every 5th label

    plt.tight_layout()
    # Show the plot
    plt.show()
    #plt.savefig(f"{Output_dir}/contact_map.png", dpi = 300)
    return residues1, residues2











def dssp_analysis(chain_id, trajectory_file, topology_file, residues_names, Output_dir, jump, name_file):
    """
    Analyze the secondary structure of a protein chain over time using DSSP and visualize the dominant structure per residue.
    
    Parameters:
    chain_id (int): The chain ID of the protein to analyze.
    trajectory_file (str): Path to the trajectory file.
    topology_file (str): Path to the topology file.
    residues_names (list of str): List of residue names for labeling the x-axis.
    
    Returns:
    None: Displays a bar plot of the most frequent secondary structure per residue.
    """
    # Load trajectory and select only protein atoms from the given chain
    traj = md.load(trajectory_file, top=topology_file)
    protein_atoms_indices = traj.topology.select(f'chainid {chain_id}')
    traj_protein_slice = traj.atom_slice(protein_atoms_indices) 
    
    # Compute secondary structure using DSSP
    dssp = md.compute_dssp(traj_protein_slice)
    
    # Classify secondary structures
    helix_codes = {'H', 'G', 'I'}  # Helices
    beta_codes = {'E', 'B'}        # Beta sheets
    coil_codes = {'T', 'S', 'C'}   # Random coil
    
    num_residues = dssp.shape[1]  # Number of residues
    
    # Count frequencies of each structure per residue
    helix_freq = np.zeros(num_residues)
    beta_freq = np.zeros(num_residues)
    coil_freq = np.zeros(num_residues)
    
    for i in range(num_residues):
        res_dssp = dssp[:, i]  # Secondary structure of residue i over time
        helix_freq[i] = np.sum(np.isin(res_dssp, list(helix_codes))) / len(res_dssp)
        beta_freq[i] = np.sum(np.isin(res_dssp, list(beta_codes))) / len(res_dssp)
        coil_freq[i] = np.sum(np.isin(res_dssp, list(coil_codes))) / len(res_dssp)
    
    # Determine the most frequent conformation
    max_freq = np.maximum(np.maximum(helix_freq, beta_freq), coil_freq)
    colors = [
        'blue' if h == mf else 'orange' if b == mf else 'red'
        for h, b, mf in zip(helix_freq, beta_freq, max_freq)
    ]
    
    # Plot
    plt.figure(figsize=(16, 8))
    plt.bar(range(num_residues), max_freq * 100, color=colors)
    plt.xlabel(f'{name_file}', fontsize=20)
    plt.ylabel('Percentage of dominant conformation', fontsize=20)
    plt.xticks(ticks=np.arange(0, len(residues_names), jump), labels=residues_names[::jump], rotation=90, fontsize=10)
    plt.yticks(fontsize=12)
    plt.gca().set_xticks(np.arange(len(residues_names)))
    
    # Add legend
    legend_handles = [
        mpatches.Patch(color='blue', label='Helix'),
        mpatches.Patch(color='orange', label='Beta sheet'),
        mpatches.Patch(color='red', label='Random coil')
    ]
    plt.legend(handles=legend_handles, loc='upper right', fontsize=18)
    plt.xticks(fontsize=12)
    plt.show()
    plt.tight_layout()
    plt.savefig(f"{Output_dir}/{name_file}dssp.png", dpi = 300)

def calculate_distance_between_centers_of_mass(u, Output_dir ,labels_size: float  = 20, ticks_size: float = 16):
    """
    This function calculates and plots the distance between the centers of mass
    of two subunits over time in a molecular dynamics simulation.

    Arguments:
    u -- MDAnalysis.Universe object representing the simulation
    """
    # Select the atoms for subunit1 and subunit2 (Cα atoms)
    subunit1 = u.select_atoms("protein and (segid A or segid seg_0 or segid seg_0_Protein_chain_A or segid Protein_chain_A) and name CA")
    subunit2 = u.select_atoms("protein and (segid B or segid seg_1 or segid seg_1_Protein_chain_B or segid Protein_chain_B) and name CA")

    # Initialize an empty list to store the distances and times
    distances = []; times = []

    # Iterate over the frames of the simulation
    for ts in u.trajectory:
        # Append the current time of the frame
        times.append(ts.time)
        # Calculate the center of mass for each subunit
        com_subunit1 = subunit1.center_of_mass()
        com_subunit2 = subunit2.center_of_mass()

        # Calculate the distance between the cember,args.x_label, args.y_label ,cutoff=8, frames_to_analyze=10)
        distance = np.linalg.norm(com_subunit1 - com_subunit2)
        distances.append(distance)

    # Convert the distances list to a numpy array for easier handling
    distances = np.array(distances)
    distances = moving_average(distances)
    plt.figure(figsize=(16, 10))
    # Plot the distances over time
    plt.plot(np.array(times)/1000, distances,linewidth=3)
    plt.xlabel('Time (ns)', fontsize = labels_size)
    plt.ylabel('Distance between Centers of Mass (Å)', fontsize = labels_size)
    #plt.title('Distance between Centers of Mass of Subunits Over Time')
    plt.yticks(fontsize=ticks_size);plt.xticks(fontsize=ticks_size)
    plt.show()
    plt.savefig(f"{Output_dir}/distance.png", dpi = 300)






def calculate_number_of_contacts(u, Output_dir, cutoff=5.0,labels_size: float  = 20, ticks_size: float = 16):
    """
    This function calculates and plots the number of contacts between two subunits
    over time in a molecular dynamics simulation.

    Arguments:
    u -- MDAnalysis.Universe object representing the simulation
    cutoff -- The distance (in Å) within which atoms are considered in contact (default is 5.0 Å)
    """
    # Select the atoms for subunit1 and subunit2 (Cα atoms)
    subunit1 = u.select_atoms("protein and (segid A or segid seg_0 or segid seg_0_Protein_chain_A or segid Protein_chain_A) and name CA")
    subunit2 = u.select_atoms("protein and (segid B or segid seg_1 or segid seg_1_Protein_chain_B or segid Protein_chain_B) and name CA")

    # Initialize lists to store the times and number of contacts
    times = []
    num_contacts = []

    # Iterate over the frames of the simulation
    for ts in u.trajectory:
        # Append the current time of the frame (convert from ps to ns)
        times.append(ts.time / 1000.0)  # Convert from ps to ns

        # Calculate the pairwise distances between atoms in subunit1 and subunit2
        contact_count = 0
        for atom1 in subunit1:
            for atom2 in subunit2:
                # Calculate the distance between the atoms
                distance = np.linalg.norm(atom1.position - atom2.position)
                if distance < cutoff:
                    contact_count += 1  # Count this as a contact if within the cutoff

        # Append the number of contacts for this frame
        num_contacts.append(contact_count)

    # Convert the lists to numpy arrays for easier handling
    times = np.array(times)
    num_contacts = np.array(num_contacts)
    num_contacts = moving_average(num_contacts)
    # Plot the number of contacts over time
    plt.figure(figsize=(16, 10))
    plt.plot(times, num_contacts,linewidth=3)
    plt.xlabel('Time (ns)', fontsize = labels_size)  # Time in nanoseconds
    plt.ylabel('Number of Contacts', fontsize = labels_size)
    plt.yticks(fontsize=ticks_size); plt.xticks(fontsize=ticks_size)
    #plt.title('Number of Contacts between Subunits Over Time', fontsize = 20)
    plt.show()
    plt.savefig(f"{Output_dir}/contacts.png", dpi = 300)

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

    print(f"Saved: {output_path}")

def plot_rmsd_rmsf(u,output_dir,residues_1_names,residues_2_names ,label_1, label_2):
    """
    Calculate and plot RMSD and RMSF for each chain in the system.

    Parameters:
    u (mda.Universe): Loaded MDAnalysis universe object containing structure and trajectory.
    """
    # RMSD calculation for each chain
    rmsd_chain = {}
      # Reference selection, usually protein backbone
    segids = set(residue.segid for residue in u.residues)
    plt.figure(figsize=(12, 6))

    for segid in segids:
        if segid == 'A': #para solo la proteina
            chain_atoms = u.select_atoms(f"segid {segid} and protein")
            print('CHAINSSSS',segid)
            # Reference atoms (positions from the first frame)
            ref_atoms = u.select_atoms(f"segid {segid} and protein")
            
            # Create the RMSD object
            rmsd_ = rms.RMSD(chain_atoms, ref_atoms)  # This will compute RMSD for each frame with respect to the first frame
            rmsd_.run()  # Run the RMSD analysis
            rmsd = rmsd_.rmsd.T
            # Plot RMSD for this chain
            #plt.plot(np.array(rmsd[1])/1000, rmsd[2], label=f"Chain {segid} RMSD")
            plt.plot(np.array(rmsd[1])/1000, rmsd[2], label=f"Protein RMSD")
    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (Å)")
    plt.legend()
    output_path = f"{output_dir}/rmsd.png"
    #plt.yticks(ticks=np.arange(0, len(residues_names), 20), labels=residues_names[::20], rotation=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()  # Close figure to avoid overlap


    
    # RMSF calculation for each chain
    count = 0
    for segid in segids:
        print('segid',segid,segids)
        try:
            try:
                plt.figure(figsize=(12, 6))
                # Select atoms for the current chain
                chain_atoms = u.select_atoms(f"segid {segid} and name CA")
                print('Trying RMSF',chain_atoms)
                print(len(chain_atoms),len(residues_1_names),len(residues_2_names))
                rmsf_values = []
                for ts in u.trajectory:
                    # Calculate the deviation for each atom
                    deviation = np.linalg.norm(chain_atoms.positions - chain_atoms.positions.mean(axis=0), axis=1)
                    rmsf_values.append(deviation)
        
                rmsf_values = np.array(rmsf_values)
        
              # Calculate the average RMSF for each atom in the chain
                avg_rmsf = rmsf_values.mean(axis=0)
        
                # Plot RMSF for the chain
                if count == 0:
                    residues_names = residues_1_names
                elif count == 1:
                    residues_names = residues_2_names
                count += 1
                print(len(residues_names),len(chain_atoms))
                plt.plot(residues_names, avg_rmsf, label=f"Chain {segid} RMSF")
                plt.xlabel("Atom Index")
                plt.ylabel("RMSF (Å)")
                plt.xticks(ticks=np.arange(0, len(residues_names), 8), labels=residues_names[::8], rotation=90, fontsize=12)
                output_path = f"{output_dir}/rmsf{segid}.png"
                #plt.yticks(ticks=np.arange(0, len(residues_names), 20), labels=residues_names[::20], rotation=0, fontsize=14)
                plt.tight_layout()
                plt.savefig(output_path, dpi=300)
                plt.close()  # Close figure to avoid overlap
            except:
                plt.figure(figsize=(12, 6))
                # Select atoms for the current chain
                chain_atoms = u.select_atoms(f"segid {segid} and name CA")
                print('Trying RMSF')
                print(len(chain_atoms),len(residues_1_names),len(residues_2_names))
                rmsf_values = []
                for ts in u.trajectory:
                    # Calculate the deviation for each atom
                    deviation = np.linalg.norm(chain_atoms.positions - chain_atoms.positions.mean(axis=0), axis=1)
                    rmsf_values.append(deviation)

                rmsf_values = np.array(rmsf_values)

                # Calculate the average RMSF for each atom in the chain
                avg_rmsf = rmsf_values.mean(axis=0)

                # Plot RMSF for the chain
                if count == 1:
                    residues_names = residues_1_names
                elif count == 0:
                    residues_names = residues_2_names
                count += 1
                print(len(residues_names),len(chain_atoms))
                plt.plot(residues_names, avg_rmsf, label=f"Chain {segid} RMSF")
                plt.xlabel("Atom Index")
                plt.ylabel("RMSF (Å)")
                plt.xticks(ticks=np.arange(0, len(residues_names), 8), labels=residues_names[::8], rotation=90, fontsize=12)
                output_path = f"{output_dir}/rmsf{segid}_chains.png"
                #plt.yticks(ticks=np.arange(0, len(residues_names), 20), labels=residues_names[::20], rotation=0, fontsize=14)
                plt.tight_layout()
                #plt.savefig(output_path, dpi=300)
                plt.close()  # Close figure to avoid overlap
        except:
            print('next')



def rmsf_2(u,residues_1_names,residues_2_names ,label_1, label_2,output_dir):
    # Cargar el universo

    # Lista de cadenas (segids)
    segids = set(residue.segid for residue in u.residues)
    # Extraer todos los segids presentes en el sistema
    segids = sorted(set(residue.segid for residue in u.residues))

    # Crear el mapa de nombres de residuos por segid
    residue_name_map = {}

    for segid in segids:
        # Seleccionar los residuos únicos del segid
        residues = u.select_atoms(f"segid {segid} and name CA").residues
        # Crear una lista de etiquetas tipo "ARG45"
        residue_names = [f"{res.resname}{i+1}" for i, res in enumerate(residues)]#[f"{res.resname}{res.resid}" for res in residues] 
        residue_name_map[segid] = residue_names
        # Alinear trayectoria completa al primer frame usando todos los CA
        aligner = align.AlignTraj(u, u, select="name CA", in_memory=True)
        aligner.run()

    # Calcular RMSF para cada cadena
    for segid in segids:
        print(f"Procesando cadena {segid}")

        # Seleccionar los átomos Cα de la cadena
        chain_atoms = u.select_atoms(f"segid {segid} and name CA")
        if len(chain_atoms) == 0:
            print(f"No se encontraron CA en segid {segid}")
            continue

        # Guardar posiciones por frame
        positions = []
        for ts in u.trajectory:
            positions.append(chain_atoms.positions.copy())
        positions = np.array(positions)  # (n_frames, n_atoms, 3)

        # Calcular RMSF
        mean_pos = positions.mean(axis=0)
        rmsf = np.sqrt(((positions - mean_pos) ** 2).sum(axis=2).mean(axis=0))

        # Recuperar los nombres de residuos desde el diccionario
        residues_names = residue_name_map.get(segid, [f"{r.resname}{r.resid}" for r in chain_atoms.residues])

        if len(residues_names) != len(rmsf):
            print(f"Longitudes incompatibles: {len(residues_names)} nombres vs {len(rmsf)} RMSF")
            continue

        # Graficar
        plt.figure(figsize=(12, 6))
        plt.plot(residues_names, rmsf, label=f"Chain {segid} RMSF")
        plt.xlabel("Residue")
        plt.ylabel("RMSF (Å)")
        plt.title(f"RMSF - Cadena {segid}")
        plt.xticks(ticks=np.arange(0, len(residues_names), 8),
                   labels=residues_names[::8], rotation=90, fontsize=10)
        plt.tight_layout()

        output_path = f"{output_dir}/rmsf_2{segid}.png"
        plt.savefig(output_path, dpi=300)
        plt.close()
        # Guardar RMSF en un archivo .txt
        txt_output_path = f"{output_dir}/rmsf_22{segid}.txt"
        with open(txt_output_path, "w") as f:
            f.write("Residue\tRMSF(Å)\n")
            for name, value in zip(residues_names, rmsf):
                f.write(f"{name}\t{value:.4f}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter specific chains from a PDB file")
    parser.add_argument("-pdb", required=True, help="Path to the input PDB file")
    parser.add_argument("-xtc", required=True, help="Trajectory file")
    parser.add_argument("-r1_number", default=0, help="Int number to start to count residue 1")
    parser.add_argument("-r2_number", default=0, help="Int number to start to count residue 2")
    parser.add_argument("-x_label", default='Residue 1', help="Name residue 1")
    parser.add_argument("-y_label", default='Residue 2', help="Name residue 2")
    parser.add_argument("-output", default='.', help="Output_dir")
    args = parser.parse_args()

    # Load trajectory once
    u = load_trajectory(args.pdb, args.xtc)


    # Perform analysis on the loaded trajectory
   
    residues_1_names, residues_2_names = contact_map(u, args.r1_number, args.r2_number,args.x_label, args.y_label, args.output, cutoff=8, frames_to_analyze=10)
    dssp_analysis(0, args.xtc, args.pdb,residues_1_names, args.output, jump=4, name_file ='Protein')
    if residues_2_names:
        dssp_analysis(1, args.xtc, args.pdb,residues_2_names, args.output, jump=1, name_file = 'Protein 2')
        dssp_chain_analysis(args.xtc, args.pdb, args.output, 'dsspchains', 1,residues_2_names)
        calculate_number_of_contacts(u, args.output, cutoff=8.0)
        calculate_distance_between_centers_of_mass(u, args.output)
    else:
        pass
    dssp_chain_analysis(args.xtc, args.pdb, args.output, 'dsspchains', 0,residues_1_names)
    plot_rmsd_rmsf(u, args.output,residues_1_names,residues_2_names, args.x_label, args.y_label)
    rmsf_2(u,residues_1_names,residues_2_names ,args.x_label, args.y_label,args.output)
