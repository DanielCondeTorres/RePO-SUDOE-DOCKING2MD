# RePO-SUDOE-DOCKING2MD

This repository automates the preparation and simulation process of protein–ligand systems from *docking* results, using **GROMACS** and auxiliary tools such as **ACPYPE** and **MDAnalysis**.

---

## 🧩 Installation

### 1️⃣ Create the Conda environment
Use the provided `.yml` file (e.g., `acpype.yml`) to create the environment on Linux:

```
conda env create -f file_name.yml
```

### 2️⃣ Activate the environment

```
conda activate acype
```

---

## 🚀 Usage

Run the complete workflow using `make`:

Makefile allows the user to just running a unique file and execute the whole code, or specific parts of it


- root: Input folder where the script will read the files                = ../Input_files/6QH4_ABEMACICLIB
- output_dir: Output folder created automatically and customized         = 6QH4_Abemaciclib_Out
- receptor_pdb: PDB file of the Protein (downloaded from the web)        = 6QH4.pdb
- ligand_name: Name of the ligand docked.                                = Abemaciclib_94NAMRd
- bfactor: Value asigned to the receptor.pdb created for VR visualization = 6QH4_Abemaciclib_94NAMRd_out.pdbqt
- resultados_vina_pdqt: Specific file coming from the RepoSUDOE Web      = 6QH4_Abemaciclib_94NAMRd_out.pdbqt.sdf
- ff: Forcefield to choose                                               = amber99sb-ildn
- water: Water model to choose                                           = tip3p


### Example Usage
```
make all \
    root="$root" \
    output_dir="$output_dir" \
    receptor_pdb="$receptor_pdb" \
    ligand_name="$ligand_name" \
    bfactor="$bfactor" \
    resultados_vina_pdqt="$resultados_vina_pdqt" \
    ff="$ff" \
    water="$water"
```


---

## ⚙️ Available Options

### Water models (`water`)
- `spc`  
- `spce`  
- `tip3p`  
- `tip4p`  
- `tip5p`  
- `tips3p`

### Force fields (`ff`)
- `charmm36-jul2022`  
- `amber99sb-ildn`  
- `amber96`  
- `charmm27`  
- `gromos54a7`  
- `gromos53a6`  
- `oplsaa`

---

## 🧰 Required Dependencies

Make sure the following packages are installed within your Conda environment (They should be):

- **GROMACS**
- **MDAnalysis**
- **mdtraj**
- **networkx**
- **prolif**

---

## 📜 Running on a Cluster (SLURM) / Interactive 

Edit the `lanzar.sh` file with your parameters:

>[!NOTE]
**Just needed to edit the $INPUT_FOLDER variable, putting the name of the folder of yours having the input files**

Then simply run:

```
sbatch lanzar.sh
```

---


## 📂 Recommended Project Structure
```
RePO-SUDOE-DOCKING2MD/
├── acpype_environment
├── Code
    ├── Working_Area
        ├── Makefile     
        ├── lanzar.sh 
        ├── Scripts_Carlos
        └── VR
    ├── Input_files
    └── Output Folders when created 
│       ...
└── README.md
```

---

## 🧪 Credits

Developed under the **SUDOE** program for workflow automation in molecular dynamics.  
Maintained by the research team **RePO-SUDOE-DOCKING2MD**.

### Acknowledgments
To Cesga, for allowing us the use of their facilities to hold this seminar.

This work has received financial support from the Spanish Agencia Estatal de Investigación (AEI) and the European Regional Development Fund - ERDF (RTI2018-098795-A-I00, PID2019-111327GB-I00, PDC2022-133402-I00 and RYC-2016-20335), by the (FPU22/00636) and the RePo-SUDOE project.

### Team
Daniel Conde-Torres, Alejandro Seco, Carlos Lozano, António Costa, Hugo A. L. Filipe, Ángel Piñeiro y Rebeca García Fandiño

### Contact
danielconde.torres@usc.es

alejandro.seco.gonzalez@usc.es

---
