# Workflow for validation of rings' conformation in ligands

## Description
* TODO:
## Getting Started
* TODO:
### Dependencies
* TODO:
### Installing
* TODO:
### Executing program
* TODO: multithreading 
* Step 1. Run PrepareDataset.py. The script gets the data for all three types of the rings at once.
It accepts two arguments:
        --pdb_local (-d): Path to the local PDB
        --pq_cmd (-p): Path to the Pattern Query .exe file

PQ Command Line version (last version 1.1.23.12.27) is included in this project.
Source: https://webchem.ncbr.muni.cz/Platform/PatternQuery

```py
python3 PrepareDataset.py -d path/to/local/pdb -p PatternQuery_1.1.23.12.27b/WebChemistry.Queries.Service.exe
```
* Step 2. Run FilterDataset.py. It accepts one argument:
        --ring (-r): Ring type (cylohexane | cyclopentane | benzene)

```py
python3 FilterDataset.py -r <ring>
```
* Step 3. Run CreateTemplates. It accepts one argument:
        --ring (-r): Ring type (cylohexane | cyclopentane | benzene)

```py
python3 CreateTemplates.py -r <ring>
```
* Step 4. Run ConfComparer.py
```py
python3 ConfComparer.py -t validation_data/<ring>/templates -i validation_data/<ring>/filtered_ligands -o validation_data/<ring>/output -e SB_batch/SiteBinderCMD.exe -d 1.0
```
* Step 5. Run analysis of electron density coverage:
```
cd electron_density_coverage_analysis
conda env create -f environment.yaml
conda activate ed-coverage-analysis
python main.py tests/example_input/validation_data -m
```
* Step 6. Run find_conformation.py.
```
cd ../validation_data/<ring>/output

```
## Help

* TODO:
```
command to run if program contains helper info
```

## Authors

Contributors names and contact info

* TODO:

## Version History
* TODO:

## License
* TODO:
This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments

* TODO: