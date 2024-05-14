# Workflow for validation of rings' conformation in ligands

## Description
This repository contains software workflow from the article by Bučeková *et al*. If you found this workflow, or any of its parts, to be useful, please cite the following article:

>BUČEKOVÁ Gabriela, Viktoriia DOSHCHENKO, Aliaksei CHARESHNEU, Jana PORUBSKÁ, Michal PAJTINKA, Michal OLEKSIK, Tomáš SVOBODA, Vladimír HORSKÝ, and Radka SVOBODOVÁ. *In preparation*. 2024.

## Prerequisites
- Python 3.10 available in the PATH environment
  - a binary installer is available from [here](https://www.python.org/downloads/release/python-31011/)
- .NET 4.0 or newer (Windows only)
- Mono 6.12.0.200 or newer (non-Windows operating systems only)
  - installation instructions are available [here](https://www.mono-project.com/download/stable)
- git (or a compatible alternative) (only if you wish to clone the workflow's repository)
- 600 GB of free space for input data

## Getting Started

1. **Clone the Project Repository (Recommended)**:
   - Open a terminal or command prompt.
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
   - Open a terminal or command prompt inside the project directory (rings-conformation-validation).

2. **Set up Virtual Environment (Optional but Recommended)**:
    - Create a virtual environment named `.venv`:
     ```
     python -m venv .venv
     ```
   - Activate the virtual environment:
     - **Windows Command Prompt**:
       ```
       .venv\Scripts\activate
       ```
     - **Windows Git Bash / macOS / Linux**:
       ```
       source .venv/bin/activate
       ```

3. **Install Required Packages**:
   - Make sure your virtual environment is activated.
   - Run:
     ```
     pip install requests gemmi==0.6.5 pandas==2.2.2 xlsxwriter requests==2.31.0

     ```

## Executing the workflow

1. Create or choose a directory where the input data will be downloaded (e.g. **user_input_dir**).
2. Create or choose a directory where the output of the workflow will be stored (e.g. **user_output_dir**).
3. Make sure you are in the project root directory **(rings-conformation-validation)**.
4. Execute the workflow:
	- Windows:
		```
		run_workflow.bat user_input_dir user_output_dir
		```
	- Linux:
		```
		bash run_workflow.sh user_input_dir user_output_dir
		```

## Using a small dataset to test the workflow

If you want to test the workflow on a small dataset, execute the workflow in the same way, as described earlier, but add the -testing switch in step 4. With the switch, the workflow execution commands are:

- Windows
	```
	run_workflow.bat -testing user_input_dir user_output_dir
	```
- Linux:
	```
	bash run_workflow.sh -testing user_input_dir user_output_dir
	```

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/sb-ncbr/rings-conformation-validation/blob/main/LICENSE) file for details.
