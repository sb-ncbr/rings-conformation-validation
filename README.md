# Workflow for validation of rings' conformation in ligands

## Description
This repository contains software workflow from the article by Bučeková et al. If you found this workflow, or any of its parts, to be useful, please cite the following article:

>BUČEKOVÁ Gabriela, Viktoriia DOSHCHENKO, Aliaksei CHARESHNEU, Jana PORUBSKÁ, Michal PAJTINKA, Michal OLEKSIK, Vladimír HORSKÝ, and Radka SVOBODOVÁ. *In preparation*. 2024.

## Prerequisites
- Python 3.10 available in the PATH environment
  - a binary installer is available from [here](https://www.python.org/downloads/release/python-31011/)
- .NET 4.0 or newer (Windows only)
- Mono 6.12.0.200 or newer (non-Windows operating systems only)
  - installation instructions are available [here](https://www.mono-project.com/download/stable)
- git (or a compatible alternative) (only if you wish to clone the workflow's repository)

## Getting Started

1. **Clone the Project Repository (Recommended)**:
   - Open a terminal or command prompt.
   - Navigate to the directory where you want to clone the project.
   - Run:
     ```
     git clone git@github.com:sb-ncbr/rings-conformation-validation.git
     ```
   - Navigate into the cloned project directory:
     ```
     cd rings-conformation-validation
     ```
   
   **Alternative: Download the workflow repository as a ZIP archive**:
   - If you're not familiar with Git or GitHub, you can download the whole workflow repository as a ZIP archive from [here](https://github.com/sb-ncbr/rings-conformation-validation/archive/refs/heads/main.zip).
   - Once downloaded, extract the ZIP archive to a directory on your computer.
   - Open a terminal or command prompt inside the project folder (rings-conformation-validation).

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
     pip install requests gemmi==0.6.5 pandas==2.2.2 xlsxwriter

     ```

## Executing program

- Create a folder where data will be downloaded (it is recommended to store the data outside the project)
- We will refer to the folder you have created as **"user_input_dir"** in further steps.
- Make sure you are in the project root directory **(rings-conformation-validation)**
- We will refer to the output folder as **"user_output_dir"**. Do not create it yourself.

### Recommended:
**Windows**
- Usage: run_workflow.bat path/to/"user_input_dir" path/to/"user_output_dir"

**Unix**
- Usage: bash run_workflow.sh path/to/"user_input_dir" path/to/"user_output_dir"

### Alternative:

0. **Run DownloadData.py**:
    - Use File ID for the "input_data" directory in OneData as a positional argument
    -  [LINK TO ONEDATA](https://doi.org/10.58074/hy79-qc22)
    - TODO: change FILE ID for real dataset in a code snippet below
        ```
        python DownloadData.py -d path/to/"user_input_dir" 00000000007ECDBE736861726547756964236133623365636535626131323133323532303238353237323438623439316133636864663563233432653234313133616330396634323834666630656235313763306539656131636865613232233630363962316339633839646164616332666562373139383633633437653639636862623462
        ```

1. **Run PrepareDataset.py**:
      - It prepares all datasets for three ring types at once.
        ```
        python PrepareDataset.py -i path/to/"user_input_dir"/"input_data" -o path/to/"user_output_dir"
        ```
2. **Run FilterDataset.py**:
    - Run this script once for each ring type.
    - Ring types: (cylohexane | cyclopentane | benzene)

        ```
        python FilterDataset.py -r <ring> -i path/to/"user_input_dir"/"input_data" -o path/to/"user_output_dir"
        ```
3. **Run CreateTemplates**:
    - Run this script once for each ring type.
    - Ring types: (cylohexane | cyclopentane )

        ```
        python CreateTemplates.py -r <ring> -i path/to/"user_input_dir"/"input_data" -o path/to/"user_output_dir"
        ```
4. **Run ConfComparer.py**:
    - Run this script once for each ring type: (cylohexane | cyclopentane | benzene)
    - For UNIX systems there is an extra argument: "-c mono" (cross-platform runner). Mono is required.
        ```
        python ConfComparer.py -t <output_dir>/validation_data/<ring>/templates -i <output_dir>/validation_data/<ring>/filtered_ligands -o <output_dir>/validation_data/<ring>/output
        ```

 5. **Run analysis of electron density coverage**:
   
        ```
        python electron_density_coverage_analysis/main.py path/to/"user_output_dir" path/to/"user_input_dir"/"input_data"/ccp4
        ```
6. **Run RingAnalysisResult**:
    - Run this script once for each ring type: (cylohexane | cyclopentane | benzene)

        ```
        python RingAnalysisiResult.py -r <ring> -i path/to/"user_input_dir"/"input_data" -o path/to/"user_output_dir"
        ```

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/sb-ncbr/rings-conformation-validation/blob/main/LICENSE) file for details.
