# Workflow for validation of rings' conformation in ligands

## Description
* TODO: update the link to the article!!!*
* This repository contains the software which is a supplementary material for [this article](https://www.overleaf.com/project/65d8603c4a4af517eb1fd627).
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
   
   **Alternative: Downloading the Project Archive**:
   - If you're not familiar with Git or GitHub, you can download the project as a ZIP archive.
     - [Download](https://github.com/sb-ncbr/rings-conformation-validation/archive/refs/heads/main.zip)
     - Once downloaded, extract the ZIP archive to a directory on your computer.
     - Open a terminal or command prompt inside the project folder (rings-conformation-validation).

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
     
### Download the PDB copy
   * Download the pdb dataset archive from [HERE](https://doi.org/10.58074/hy79-qc22) and extract it into the folder of your choice on your computer

### Download CCP4 files
   * Change the directory to "electron_density_coverage_analysis"
   * Create a folder "ccp4" inside this directory
   * Download the ccp4 archive from [HERE](https://doi.org/10.58074/hy79-qc22) and extract it into "ccp4" folder you have just created
    
### Download Chemical Component Dictionary
   * Download components.cif.gz into the project root folder "rings-conformation-validation". Do not extract the file from the archive.


### Executing program

1. **Run PrepareDataset.py**:
    - Navigate into project root folder "rings-conformation-validation".
    - Before proceeding, make sure you have downloaded the copy of the PDB into your local machine. [Jump to instructions](#download-the-pdb-copy)
    - Make sure you have components.cif.gz in your current working directory [Jump to instructions](#download-chemical-component-dictionary)
    - [PQ Command Line version](https://webchem.ncbr.muni.cz/Platform/PatternQuery) (last version 1.1.23.12.27) is included in this project.
      
        ```
        python PrepareDataset.py -d path/to/local/pdb
        ```
2. **Run FilterDataset.py**:
    - It accepts one argument --ring (-r): Ring type (cylohexane | cyclopentane | benzene)

        ```
        python FilterDataset.py -r <ring>
        ```
3. **Run CreateTemplates**:
    - It accepts one argument --ring (-r): Ring type (cylohexane | cyclopentane | benzene)

        ```
        python CreateTemplates.py -r <ring>
        ```
4. **Run ConfComparer.py**:
    - For benzene ring run the script twice, fisrt time normally, the second time use templates from cyclohexane
    - TODO: correct text
    
        ```
        python ConfComparer.py -t validation_data/<ring>/templates -i validation_data/<ring>/filtered_ligands -o validation_data/<ring>/output
        ```
 5. **Run analysis of electron density coverage**:
    - Change the directory to <electron_density_coverage_analysis>
    
        ```py
        cd electron_density_coverage_analysis
        ```
    - Before proceeding, make sure you have downloaded the CCP4 files into "ccp4" folder. [Jump to instructions](#download-ccp4-files)
   
        ```
        python main.py validation_data
        ```
6. **Run RingAnalysisResult**:
    - It accepts one argument --ring (-r): Ring type (cylohexane | cyclopentane | benzene)

        ```
        python RingAnalysisiResult.py -r <ring>
        ```

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/sb-ncbr/rings-conformation-validation/blob/main/LICENSE) file for details.
