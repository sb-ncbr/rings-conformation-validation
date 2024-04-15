# Workflow for validation of rings' conformation in ligands

## Description
* TODO:
## Getting Started

1. **Clone the Project Repository**:
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

2. **Install Python 3.10**: Download and install Python 3.10 from the official website: [Python Downloads](https://www.python.org/downloads/)
   
   **Windows and macOS**: Ensure to check the option to add Python to the PATH during installation.

3. **Set up Virtual Environment (Optional but Recommended)**:
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

4. **Install Required Packages**:
   - Make sure your virtual environment is activated.
   - Run:
     ```
     pip install gemmi==0.6.5 pandas==2.2.1

     ```

5. **Install External Tools**:
   - **Non-Windows Operating System (Linux and macOS)**:
     - Make sure you have installed Mono. You can download and install it from [Mono Project](http://mono-project.com).
     - Current LST version is 6.12.0 
   
   - **Windows Operating System**:
     - Make sure you have the .NET environment installed (version 4.0 and above). By default, they are on all up-to-date Windows operating systems.

### Executing program

* Step 1. Run PrepareDataset.py. The script gets the data for all three types of the rings at once.
It accepts two arguments:
        --pdb_local (-d): Path to the local PDB
        --pq_cmd (-p): Path to the Pattern Query .exe file

PQ Command Line version (last version 1.1.23.12.27) is included in this project.
Source: https://webchem.ncbr.muni.cz/Platform/PatternQuery

```py
python PrepareDataset.py -d path/to/local/pdb -p PatternQuery_1.1.23.12.27b/WebChemistry.Queries.Service.exe
```
* Step 2. Run FilterDataset.py. It accepts one argument:
        --ring (-r): Ring type (cylohexane | cyclopentane | benzene)

```py
python FilterDataset.py -r <ring>
```
* Step 3. Run CreateTemplates. It accepts one argument:
        --ring (-r): Ring type (cylohexane | cyclopentane | benzene)

```py
python CreateTemplates.py -r <ring>
```
* Step 4. Run ConfComparer.py
```py
python ConfComparer.py -t validation_data/<ring>/templates -i validation_data/<ring>/filtered_ligands -o validation_data/<ring>/output -e SB_batch/SiteBinderCMD.exe -d 1.0
```
* Step 5. Run analysis of electron density coverage:
```
cd electron_density_coverage_analysis
python main.py validation_data -m
```
* Step 6. Run RingAnalysisResult. It accepts one argument:
        --ring (-r): Ring type (cylohexane | cyclopentane | benzene)
* TODO: correct paths
```
cd ../validation_data/<ring>/output
python RingAnalysisiResult.py -r <ring>

```
## Help

* TODO:
```
command to run if program contains helper info
```

## Authors

Contributors names and contact info

* TODO:


## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

* TODO: