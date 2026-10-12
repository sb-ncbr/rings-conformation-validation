# Workflow for validation of rings' conformation in non-standard residues

## Description
The original software workflow is from the article by Bučeková *et al.* *Analysis of cyclopentane, cyclohexane, and benzene conformations in ligands for PDB x-ray structures using the Hill-Reilly approach* that has been submitted to the Journal of Cheminformatics.

Current release contains an updated software workflow from the article by Doshchenko *et al.* *RingPDB: a web application for exploration of ring conformations in non-standard residues from the Protein Data Bank* that has been submitted to the Journal of Cheminformatics.


## Prerequisites
- Python 3.12.3 available in the PATH environment
	- a binary installer is available from [here](https://www.python.org/downloads/release/python-3123/)
- git (or a compatible alternative)
- A local copy of coordinate files in PDBx/mmCIF Format retrieved from the [PDB Beta Archive](https://www.wwpdb.org/ftp/pdb-beta-ftp-sites). It is important to maintain the archive directory hierarchy.
- A locally available copy of the [Chemical Component Dictionary (CCD)](https://doi.org/10.1093/bioinformatics/btu789).
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
	- Create and activate virtual environment named `.venv`:

		**macOS / Linux:**
		```
		python3 -m venv .venv
		source .venv/bin/activate
		```
		**Windows:**
		```
		python -m venv .venv
		.venv\Scripts\activate
		```

3. **Install Required Packages**:
	- Make sure your virtual environment is activated.
	- Run:
		```
		python -m pip install -r requirements.txt
		```

## Executing the workflow

1. Create `config.yaml` file and copy the content of `config_example.yaml` into your file. Change the paths according to the description in the yaml file.
4. Execute the workflow:
	```
	python run.py
	```

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/sb-ncbr/rings-conformation-validation/blob/main/LICENSE) file for details.
