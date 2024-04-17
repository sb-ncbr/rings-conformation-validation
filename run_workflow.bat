@echo off

set "LOCAL_PDB_PATH=%~1"
set "OUTPUT_FOLDER=%~2"

set "CYCLOPENTANE=cyclopentane"
set "CYCLOHEXANE=cyclohexane"
set "BENZENE=benzene"
set list=cyclopentane cyclohexane benzene


REM python PrepareDataset.py -d "%LOCAL_PDB_PATH%" -o "%OUTPUT_FOLDER%"

for %%r in (%CYCLOPENTANE% %CYCLOHEXANE% %BENZENE%) do (
    REM python FilterDataset.py -r "%%~r" -o "%OUTPUT_FOLDER%"

    python CreateTemplates.py -r "%%~r" -o "%OUTPUT_FOLDER%"

    set "TEMPLATES_PATH=%OUTPUT_FOLDER%\validation_data\"%%~r"\templates"
    set "FILTERED_LIGANDS_PATH=%OUTPUT_FOLDER%\validation_data\%RING%\filtered_ligands"
    set "OUTPUT_PATH=%OUTPUT_FOLDER%\validation_data\%RING%\output"

    python ConfComparer.py -t "%TEMPLATES_PATH%" -i "%FILTERED_LIGANDS_PATH%" -o "%OUTPUT_PATH%"
)


cd electron_density_coverage_analysis
python main.py %OUTPUT_FOLDER%\validation_data

REM Step 6: Run RingAnalysisResult.py
REM python RingAnalysisResult.py -r <ring>
