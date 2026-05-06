#!/bin/bash

CYCLOPENTANE="cyclopentane"
CYCLOHEXANE="cyclohexane"
BENZENE="benzene"
OXANE="oxane"
OXOLANE="oxolane"

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
# Note: $INPUT_DATA_FOLDER should be created before running this script
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

python3 FilterDataset.py -i "$INPUT_DATA_FOLDER/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "Error: FilterDataset failed with exit code $exit_code"
    exit $exit_code
fi

python3 CalculateHR.py "$OUTPUT_FOLDER"
python3 CompareHR.py "$OUTPUT_FOLDER"

# analyse electron density coverage
python3 AnalyseCoverage.py "$OUTPUT_FOLDER" "${INPUT_DATA_FOLDER}/${DATA_FOLDER}/ccp4"

# analyse and summarise results
python3 RingAnalysisResult.py -r "$CYCLOPENTANE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
python3 RingAnalysisResult.py -r "$CYCLOHEXANE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
python3 RingAnalysisResult.py -r "$BENZENE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
python3 RingAnalysisResult.py -r "$OXANE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
python3 RingAnalysisResult.py -r "$OXOLANE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"

python3 BuildWebDataset.py -i "$INPUT_DATA_FOLDER/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "Error: BuildWebDataset failed with exit code $exit_code"
    exit $exit_code
fi