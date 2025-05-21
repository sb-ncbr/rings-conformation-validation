# Workflow for validation of rings' conformation in ligands

## Description
This repository contains the software workflow from the article by Bučeková *et al.* *Analysis of cyclopentane, cyclohexane and benzene conformations in ligands for PDB x-ray structures* that has been submitted to the Journal of Cheminformatics.

## Prerequisites
- Python 3.12.3 available in the PATH environment
	- a binary installer is available from [here](https://www.python.org/downloads/release/python-3123/
- Mono 6.12.0.200 or newer (non-Windows operating systems only)
	- installation instructions are available [here](https://www.mono-project.com/download/stable
- git (or a compatible alternative) (only if you wish to clone the workflow's repository)
- 1 TB of free space for input and output data

## Getting Started

1. **Clone the Project Repository (Recommended)**:
	- Open a terminal or a command prompt.
	- Navigate to the directory where you want to clone the project.
	- Run:
	```
	git clone https://github.com/sb-ncbr/rings-conformation-validation.git
	```
	- Navigate into the cloned project directory:
	```
	cd rings-conformation-validation
	```
	
	**Alternative: Download the workflow repository as a ZIP archive**:
	- If you're not familiar with Git or GitHub, you can download the whole workflow repository as a ZIP archive from [here](https://github.com/sb-ncbr/rings-conformation-validation/archive/refs/heads/main.zip).
	- Once downloaded, extract the ZIP archive to a directory on your computer.
	- Open a terminal or a command prompt inside the project directory (rings-conformation-validation).

2. **Set up Virtual Environment (Optional but Recommended)**:
	- Create a virtual environment named `.venv`:
	```
	python3 -m venv .venv
	```
	- Activate the virtual environment:
	```
	source .venv/bin/activate
	```

3. **Install Required Packages**:
	- Make sure your virtual environment is activated.
	- Run:
	```
	pip install gemmi==0.6.5 pandas==2.2.2 xlsxwriter requests==2.32.3 Scipy==1.15.2 Numpy==2.1.3 Numba==0.61.0 Biopython==1.85
	```

## Executing the workflow

1. Create or choose a directory where the input data will be downloaded (e.g. **user_input_dir**).
2. Create or choose a directory where the output of the workflow will be stored (e.g. **user_output_dir**).
3. Make sure you are in the project root directory **(rings-conformation-validation)**.
4. Execute the workflow:
	```
	bash run_workflow.sh user_input_dir user_output_dir
	```

## Using a small dataset to test the workflow

If you want to test the workflow on a small dataset, execute the workflow in the same way as described earlier, but add the -testing switch in step 4. With the switch, the workflow execution commands are:

```
bash run_workflow.sh -testing user_input_dir user_output_dir
```

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/sb-ncbr/rings-conformation-validation/blob/main/LICENSE) file for details.
