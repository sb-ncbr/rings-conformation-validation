@echo off
setlocal enabledelayedexpansion

set "CYCLOPENTANE=cyclopentane"
set "CYCLOHEXANE=cyclohexane"
set "BENZENE=benzene"
set "ONEDATA_ID=000000000052BCFF67756964233134623266393865646664313231356137333837626637373163306234633530636866363339233432653234313133616330396634323834666630656235313763306539656131636865613232"
set "DATA_FOLDER=input_data"

if "%1" == "-testing" (
    echo Testing is ON
    set "DATA_FOLDER=input_data_small"
    set "ONEDATA_ID=000000000052D6A467756964236133623365636535626131323133323532303238353237323438623439316133636864663563233432653234313133616330396634323834666630656235313763306539656131636865613232"
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
