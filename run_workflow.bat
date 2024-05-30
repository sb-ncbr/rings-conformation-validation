@echo off
setlocal enabledelayedexpansion

set "CYCLOPENTANE=cyclopentane"
set "CYCLOHEXANE=cyclohexane"
set "BENZENE=benzene"
set "ONEDATA_ID=00000000007E7B78736861726547756964236135336261663062373836396337623335393963643438393464303838393265636866656162236563663834616464323164326666613165373037633331393464326264633264636830303830236438383464313936663731323531356638653938316631313833636338363666636835313463"
set "DATA_FOLDER=input_data"

if "%1" == "-testing" (
    echo Testing is ON
    set "DATA_FOLDER=input_data_small"
    set "ONEDATA_ID=00000000007E3678736861726547756964236232663464646532663737336664626564653862386165373530373631346139636838353230236563663834616464323164326666613165373037633331393464326264633264636830303830236438383464313936663731323531356638653938316631313833636338363666636835313463"
    shift
)

set "INPUT_DATA_FOLDER=%~1"
set "OUTPUT_FOLDER=%~2"

REM $INPUT_DATA_FOLDER should be created
python DownloadData.py -j 4 -d !INPUT_DATA_FOLDER! !ONEDATA_ID!
if "!errorlevel!" NEQ "0" (
     exit /b !errorlevel!
)

python PrepareDataset.py -i !INPUT_DATA_FOLDER!\!DATA_FOLDER! -o !OUTPUT_FOLDER!

if "!errorlevel!" NEQ "0" (
    exit /b !errorlevel!
)

for %%r in (%CYCLOPENTANE%, %CYCLOHEXANE%) do (

    python FilterDataset.py -r %%~r -i !INPUT_DATA_FOLDER!\!DATA_FOLDER! -o !OUTPUT_FOLDER!
    if "!errorlevel!" NEQ "0" (
        exit /b !errorlevel!
    )

    python CreateTemplates.py -r %%~r -i !INPUT_DATA_FOLDER!\!DATA_FOLDER! -o !OUTPUT_FOLDER!
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

python FilterDataset.py -r !BENZENE! -i !INPUT_DATA_FOLDER!\!DATA_FOLDER! -o !OUTPUT_FOLDER!
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

set "CCP4=!INPUT_DATA_FOLDER!\!DATA_FOLDER!\ccp4"

python electron_density_coverage_analysis\main.py !OUTPUT_FOLDER! !CCP4!
if "!errorlevel!" NEQ "0" (
    exit /b !errorlevel!
)

python RingAnalysisResult.py -r !CYCLOPENTANE! -i !INPUT_DATA_FOLDER!\!DATA_FOLDER! -o !OUTPUT_FOLDER!

python RingAnalysisResult.py -r !CYCLOHEXANE! -i !INPUT_DATA_FOLDER!\!DATA_FOLDER! -o !OUTPUT_FOLDER!

python RingAnalysisResult.py -r !BENZENE! -i !INPUT_DATA_FOLDER!\!DATA_FOLDER! -o !OUTPUT_FOLDER!

endlocal
