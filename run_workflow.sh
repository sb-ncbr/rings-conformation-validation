#!/bin/bash

CYCLOPENTANE="cyclopentane"
CYCLOHEXANE="cyclohexane"
BENZENE="benzene"
ONEDATA_ID=00000000007EFEE3736861726547756964236235336339646566653264323965303965386534313264623166333864373064636832353434236563663834616464323164326666613165373037633331393464326264633264636830303830236437616333386562636534633830396363316166343935303363346236333962636838353261
DATA_FOLDER="input_data"

if [[ "$1" == "-testing" ]]; then
    echo "Testing is ON"
    DATA_FOLDER="input_data_small"
    ONEDATA_ID=00000000007EB4F7736861726547756964233735373966346537353231383834613064333339373764313438303731346534636831643962236563663834616464323164326666613165373037633331393464326264633264636830303830236437616333386562636534633830396363316166343935303363346236333962636838353261
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

# Download of data
# None: $INPUT_DATA_FOLDER should be created before running this script
python3 DownloadData.py -j 4 -d "$INPUT_DATA_FOLDER" $ONEDATA_ID
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "Error: DownloadData failed with exit code $exit_code"
    exit $exit_code
fi

# Prepare dataset using PatternQuery
python3 PrepareDataset.py -i "$INPUT_DATA_FOLDER/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "Error: PrepareDataset failed with exit code $exit_code"
    exit $exit_code
fi

# Identify conformation of cyclohexane cycles
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

# Identify conformation of cyclopentane cycles
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

# Identify conformation of benzene cycles
# Note: in directory QM_optimised_templates/benzene is file flat.pdb for benzene and another conformations for cyclohexane
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

# analyse electron density coverage
CCP4="${INPUT_DATA_FOLDER}/${DATA_FOLDER}/ccp4"
python3 electron_density_coverage_analysis/main.py "$OUTPUT_FOLDER" "$CCP4"

# analyse and summarise results
python3 RingAnalysisResult.py -r "$CYCLOPENTANE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
python3 RingAnalysisResult.py -r "$CYCLOHEXANE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
python3 RingAnalysisResult.py -r "$BENZENE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"