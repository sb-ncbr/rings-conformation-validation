@echo off
setlocal enabledelayedexpansion
set "INPUT_DATA_FOLDER=%~1"
set "OUTPUT_FOLDER=%~2"

set "CYCLOPENTANE=cyclopentane"
set "CYCLOHEXANE=cyclohexane"
set "BENZENE=benzene"

REM TODO: update later!
REM for testing only
set "DATA_FOLDER=input_data_small"
REM set "DATA_FOLDER=input_data"

REM $INPUT_DATA_FOLDER should be created
python DownloadData.py -d !INPUT_DATA_FOLDER! 00000000007ECDBE736861726547756964236133623365636535626131323133323532303238353237323438623439316133636864663563233432653234313133616330396634323834666630656235313763306539656131636865613232233630363962316339633839646164616332666562373139383633633437653639636862623462

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
