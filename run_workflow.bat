@echo off
setlocal enabledelayedexpansion

set "CYCLOPENTANE=cyclopentane"
set "CYCLOHEXANE=cyclohexane"
set "BENZENE=benzene"
set "ONEDATA_ID=000000000052E8FC67756964233464666266383736383334666338643632316565356135313136316239313733636866613736233536353132643962373637363232356635353333643661396662333364363936636839303432"
set "DATA_FOLDER=input_data"

if "%1" == "-testing" (
    echo Testing is ON
    set "DATA_FOLDER=input_data_small"
    set "ONEDATA_ID=000000000052713767756964236664613366353531383739356131343932336331663538386138663466383533636866616565233536353132643962373637363232356635353333643661396662333364363936636839303432"
    shift
)

set "INPUT_DATA_FOLDER=%~1"
set "OUTPUT_FOLDER=%~2"

REM $INPUT_DATA_FOLDER should be created
 python DownloadData.py -d !INPUT_DATA_FOLDER! !ONEDATA_ID!
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
