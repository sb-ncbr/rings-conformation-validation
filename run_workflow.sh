#!/bin/bash

CYCLOPENTANE="cyclopentane"
CYCLOHEXANE="cyclohexane"
BENZENE="benzene"
ONEDATA_ID=00000000007E7B78736861726547756964236135336261663062373836396337623335393963643438393464303838393265636866656162236563663834616464323164326666613165373037633331393464326264633264636830303830236438383464313936663731323531356638653938316631313833636338363666636835313463
DATA_FOLDER="input_data"

if [[ "$1" == "-testing" ]]; then
    echo "Testing is ON"
    DATA_FOLDER="input_data_small"
    ONEDATA_ID=00000000007E3678736861726547756964236232663464646532663737336664626564653862386165373530373631346139636838353230236563663834616464323164326666613165373037633331393464326264633264636830303830236438383464313936663731323531356638653938316631313833636338363666636835313463

    shift
fi


usage() {
    echo "Usage: $0 [-testing] <user_input_dir> <user_output_dir>"
    exit 1
}

if [ "$#" -ne 2 ]; then
    usage
fi

INPUT_DATA_FOLDER="$1"
OUTPUT_FOLDER="$2"

#$INPUT_DATA_FOLDER should be created
python3 DownloadData.py -j 4 -d "$INPUT_DATA_FOLDER" $ONEDATA_ID
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "Error: DownloadData failed with exit code $exit_code"
    exit $exit_code
fi

python3 PrepareDataset.py -i "$INPUT_DATA_FOLDER/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "Error: PrepareDataset failed with exit code $exit_code"
    exit $exit_code
fi




# cyclohexane
python3 FilterDataset.py -r "cyclohexane" -i "$INPUT_DATA_FOLDER/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "Error: FilterDataset failed with exit code $exit_code"
    exit $exit_code
fi
FILTERED_LIGANDS_PATH="$OUTPUT_FOLDER/validation_data/cyclohexane/filtered_ligands"
CC_OUTPUT_PATH="$OUTPUT_FOLDER/validation_data/cyclohexane/output"
mkdir "$CC_OUTPUT_PATH"
python3 SelectConformation.py "cyclohexane" "$FILTERED_LIGANDS_PATH" "$CC_OUTPUT_PATH"


# cyclopentane
python3 FilterDataset.py -r "cyclopentane" -i "$INPUT_DATA_FOLDER/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "Error: FilterDataset failed with exit code $exit_code"
    exit $exit_code
fi
FILTERED_LIGANDS_PATH="$OUTPUT_FOLDER/validation_data/cyclopentane/filtered_ligands"
CC_OUTPUT_PATH="$OUTPUT_FOLDER/validation_data/cyclopentane/output"
mkdir "$CC_OUTPUT_PATH"
python3 SelectConformation.py "cyclopentane" "$FILTERED_LIGANDS_PATH" "$CC_OUTPUT_PATH"


# benzene
# in directory QM_optimised_templates/benzene is file flat.pdb for benzene and another conformations for cyclohexane
python3 FilterDataset.py -r "benzene" -i "$INPUT_DATA_FOLDER/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "Error: FilterDataset failed with exit code $exit_code"
    exit $exit_code
fi
FILTERED_LIGANDS_PATH="$OUTPUT_FOLDER/validation_data/benzene/filtered_ligands"
CC_OUTPUT_PATH="$OUTPUT_FOLDER/validation_data/benzene/output"
mkdir "$CC_OUTPUT_PATH"
python3 SelectConformation.py "benzene" "$FILTERED_LIGANDS_PATH" "$CC_OUTPUT_PATH"


CCP4="${INPUT_DATA_FOLDER}/${DATA_FOLDER}/ccp4"

python3 electron_density_coverage_analysis/main.py "$OUTPUT_FOLDER" "$CCP4"

python3 RingAnalysisResult.py -r "$CYCLOPENTANE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
python3 RingAnalysisResult.py -r "$CYCLOHEXANE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
python3 RingAnalysisResult.py -r "$BENZENE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"