@echo off
setlocal enabledelayedexpansion
set "INPUT_DATA_FOLDER=%~1"
set "OUTPUT_FOLDER=%~2"

set "CYCLOPENTANE=cyclopentane"
set "CYCLOHEXANE=cyclohexane"
set "BENZENE=benzene"


python PrepareDataset.py -i !INPUT_DATA_FOLDER! -o !OUTPUT_FOLDER!

for %%r in (%CYCLOPENTANE%, %CYCLOHEXANE%) do (
    python FilterDataset.py -r %%~r -o !OUTPUT_FOLDER! -i !INPUT_DATA_FOLDER!
    python CreateTemplates.py -r %%~r -o !OUTPUT_FOLDER!

    set "TEMPLATES_PATH=!OUTPUT_FOLDER!\validation_data\%%~r\templates"
    set "FILTERED_LIGANDS_PATH=!OUTPUT_FOLDER!\validation_data\%%~r\filtered_ligands"
    set "OUTPUT_PATH=!OUTPUT_FOLDER!\validation_data\%%~r\output"

    python ConfComparer.py -t !TEMPLATES_PATH! -i !FILTERED_LIGANDS_PATH! -o !OUTPUT_PATH!
)

python FilterDataset.py -r !BENZENE! -o !OUTPUT_FOLDER! -i !INPUT_DATA_FOLDER!

set "TEMPLATES_PATH=!OUTPUT_FOLDER!\validation_data\!CYCLOHEXANE!\templates"
set "FILTERED_LIGANDS_PATH=!OUTPUT_FOLDER!\validation_data\!BENZENE!\filtered_ligands"
set "OUTPUT_PATH=!OUTPUT_FOLDER!\validation_data\!BENZENE!\output"

python ConfComparer.py -t !TEMPLATES_PATH! -i !FILTERED_LIGANDS_PATH! -o !OUTPUT_PATH!


REM cd electron_density_coverage_analysis
REM python main.py %OUTPUT_FOLDER%\validation_data

REM Step 6: Run RingAnalysisResult.py
REM python RingAnalysisResult.py -r <ring>

endlocal
