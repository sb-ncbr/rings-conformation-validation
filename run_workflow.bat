@echo off
setlocal enabledelayedexpansion
set "INPUT_DATA_FOLDER=%~1"
set "OUTPUT_FOLDER=%~2"

set "CYCLOPENTANE=cyclopentane"
set "CYCLOHEXANE=cyclohexane"
set "BENZENE=benzene"


python PrepareDataset.py -i !INPUT_DATA_FOLDER! -o !OUTPUT_FOLDER!

if "!errorlevel!" NEQ "0" (
    exit /b !errorlevel!
)

for %%r in (%CYCLOPENTANE%, %CYCLOHEXANE%) do (

    python FilterDataset.py -r %%~r -o !OUTPUT_FOLDER! -i !INPUT_DATA_FOLDER!
    if "!errorlevel!" NEQ "0" (
        exit /b !errorlevel!
    )

    python CreateTemplates.py -r %%~r -o !OUTPUT_FOLDER! -i !INPUT_DATA_FOLDER!
    if "!errorlevel!" NEQ "0" (
        exit /b !errorlevel!
    )

    set "TEMPLATES_PATH=!OUTPUT_FOLDER!\validation_data\%%~r\templates"
    set "FILTERED_LIGANDS_PATH=!OUTPUT_FOLDER!\validation_data\%%~r\filtered_ligands"
    set "OUTPUT_PATH=!OUTPUT_FOLDER!\validation_data\%%~r\output"

    python ConfComparer.py -t !TEMPLATES_PATH! -i !FILTERED_LIGANDS_PATH! -o !OUTPUT_PATH!
    if "!errorlevel!" NEQ "0" (
        exit /b !errorlevel!
    )
)

python FilterDataset.py -r !BENZENE! -o !OUTPUT_FOLDER! -i !INPUT_DATA_FOLDER!
if "!errorlevel!" NEQ "0" (
    exit /b !errorlevel!
)
set "TEMPLATES_PATH=!OUTPUT_FOLDER!\validation_data\!CYCLOHEXANE!\templates"
set "FILTERED_LIGANDS_PATH=!OUTPUT_FOLDER!\validation_data\!BENZENE!\filtered_ligands"
set "OUTPUT_PATH=!OUTPUT_FOLDER!\validation_data\!BENZENE!\output"

python ConfComparer.py -t !TEMPLATES_PATH! -i !FILTERED_LIGANDS_PATH! -o !OUTPUT_PATH!
if "!errorlevel!" NEQ "0" (
    exit /b !errorlevel!
)

cd electron_density_coverage_analysis
REM set "CCP4=!INPUT_DATA_FOLDER!\ccp4"
set "CCP4=.\ccp4"

python main.py !OUTPUT_FOLDER! !CCP4!
if "!errorlevel!" NEQ "0" (
    exit /b !errorlevel!
)

cd ..
python RingAnalysisResult.py -r !CYCLOPENTANE! -o !OUTPUT_FOLDER! -i !INPUT_DATA_FOLDER!

python RingAnalysisResult.py -r !CYCLOHEXANE! -o !OUTPUT_FOLDER! -i !INPUT_DATA_FOLDER!

python RingAnalysisResult.py -r !BENZENE! -o !OUTPUT_FOLDER! -i !INPUT_DATA_FOLDER!

endlocal
