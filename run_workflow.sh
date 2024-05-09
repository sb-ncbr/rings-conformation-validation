#!/bin/bash

CYCLOPENTANE="cyclopentane"
CYCLOHEXANE="cyclohexane"
BENZENE="benzene"
ONEDATA_ID=000000000052BCFF67756964233134623266393865646664313231356137333837626637373163306234633530636866363339233432653234313133616330396634323834666630656235313763306539656131636865613232
DATA_FOLDER="input_data"

if [[ "$1" == "-testing" ]]; then
    echo "Testing is ON"
    DATA_FOLDER="input_data_small"
    ONEDATA_ID=000000000052D6A467756964236133623365636535626131323133323532303238353237323438623439316133636864663563233432653234313133616330396634323834666630656235313763306539656131636865613232

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

# $INPUT_DATA_FOLDER should be created
python3 DownloadData.py -d "$INPUT_DATA_FOLDER" $ONEDATA_ID
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

for r in "$CYCLOPENTANE" "$CYCLOHEXANE"; do
    python3 FilterDataset.py -r "$r" -i "$INPUT_DATA_FOLDER/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
    exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "Error: FilterDataset failed with exit code $exit_code"
        exit $exit_code
    fi

    python3 CreateTemplates.py -r "$r" -i "$INPUT_DATA_FOLDER/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
    exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "Error: CreateTemplates failed with exit code $exit_code"
        exit $exit_code
    fi

    TEMPLATES_PATH="$OUTPUT_FOLDER/validation_data/$r/templates"
    FILTERED_LIGANDS_PATH="$OUTPUT_FOLDER/validation_data/$r/filtered_ligands"
    CC_OUTPUT_PATH="$OUTPUT_FOLDER/validation_data/$r/output"

    python3 ConfComparer.py -t "$TEMPLATES_PATH" -i "$FILTERED_LIGANDS_PATH" -o "$CC_OUTPUT_PATH" -c mono
    exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "Error: ConfComparer failed with exit code $exit_code"
        exit $exit_code
    fi
done

python3 FilterDataset.py -r "$BENZENE" -i "$INPUT_DATA_FOLDER/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "Error: FilterDataset failed with exit code $exit_code"
        exit $exit_code
    fi

# Benzenes are compared with templates for cyclohexane
TEMPLATES_PATH="$OUTPUT_FOLDER/validation_data/$CYCLOHEXANE/templates"
FILTERED_LIGANDS_PATH="$OUTPUT_FOLDER/validation_data/$BENZENE/filtered_ligands"
CC_OUTPUT_PATH="$OUTPUT_FOLDER/validation_data/$BENZENE/output"
python3 ConfComparer.py -t "$TEMPLATES_PATH" -i "$FILTERED_LIGANDS_PATH" -o "$CC_OUTPUT_PATH" -c mono
exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "Error: ConfComparer failed with exit code $exit_code"
        exit $exit_code
    fi

CCP4="${INPUT_DATA_FOLDER}/${DATA_FOLDER}/ccp4"

python3 electron_density_coverage_analysis/main.py "$OUTPUT_FOLDER" "$CCP4"

python3 RingAnalysisResult.py -r "$CYCLOPENTANE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
python3 RingAnalysisResult.py -r "$CYCLOHEXANE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"
python3 RingAnalysisResult.py -r "$BENZENE" -i "${INPUT_DATA_FOLDER}/${DATA_FOLDER}" -o "$OUTPUT_FOLDER"