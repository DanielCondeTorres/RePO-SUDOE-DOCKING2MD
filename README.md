# RePO-SUDOE-DOCKING2MD

## Install

Use the acpype .yml for Linux

```
conda env create -f nombre_del_archivo.yml
```

## Usage
```
conda activate acype
```


This workflow is automated via the `lanzar.sh` script. The process is designed to minimize user input by automatically deriving necessary simulation variables from the input files.

### 1. Configuration

To run a simulation, you only need to modify **one variable** inside the `lanzar.sh` script:

* **`INPUT_FOLDER`**: Set this to the name of the directory containing your specific input files.
* *Note:* Ensure the files inside this folder share the corresponding naming convention (root name), as the script uses these filenames to define the variables passed to the `make` command.



> **Important:** The **Forcefield** and **Water Model** parameters are currently fixed within `lanzar.sh`.

Inside lanzar.sh you only need to modify INPUT_FOLDER name (this folder must be inside de directory Input_files)
```
# =====================================
# Carpeta en Input_Files que va a usarse
# =====================================

INPUT_FOLDER="2REI_Afatinib"

# --- Definir variables automáticamente ---
ROOT="../Input_files/$INPUT_FOLDER"
eval "$("$DEF_VARS_SCRIPT" "$ROOT")"
```

### 2. Execution

Once the input folder is defined, execute the script:

```bash
./lanzar.sh

```

**Internal Logic:**

1. The script invokes `Makefile_v2` (note: the standard `Makefile` is present but ignored).
2. It parses the `INPUT_FOLDER` to configure the simulation.
3. Upon completion, it automatically runs `make clean` to remove temporary simulation artifacts.

### 3. Output Structure

The script generates a specific output directory containing the following sub-folders:

* 📂 **Docking_picture**
* 📂 **Docking_poses**
* 📂 **Dynamics_results**
* 📂 **Simulation_files**
* 📂 **VR**

**Cleanup & Logs:**

* **Temporary Files:** Any files matching the pattern `#*` are automatically deleted from the output to keep the directory clean.
* **SLURM Log:** A `slurm.out` file is generated in the `Working_Area`.
* **Process Log:** A specific log file named `Ligand_Code_Receptor.log` is generated in the output folder.
* This log tracks the execution of `Makefile_v2`, recording which targets were executed and confirming if the process finished normally or encountered errors.



---

### ⚠️ Technical Notes (For Developers)

* **`main_ligand_modified.sh`**: You may notice this script in the directory. It is an experimental version intended to test resource allocation (specifically `ntmpi` and `ntomp` combinations like 4/1, 2/2, etc.). **It is currently unused** due to potential stability issues. The active `lanzar.sh` defaults to the standard `ntmpi 1` configuration.

---


```
## Options
- Water: spc, spce, tip3p, tip4p, tip5p, tips3p
- FF:  charmm36-jul2022, amber99sb-ildn, amber96, charmm27, gromos54a7, gromos53a6, oplsaa

## Necesario:
- Gromacs
- MDanalsys
- mdtraj
- networkx
- prolif


# NOTE!!!!!:

Change in lanzar.sh 

```
make all root=../Input_files/experiment_39_docking_files/ output_dir=Prueba_final receptor_pdb=6QH4.pdb ligand_name=Palbociclib  resultados_vina_pdqt=6QH4_Palbociclib_out.pdbqt.sdf ff=amber99sb-ildn water=tip3p
```

Only is needed to use:

```
sbatch lanzar.sh
```


Necesario darle un par de vueltas al pdbfixer.py, pero es necesario para arreglar los pdb iniciales
