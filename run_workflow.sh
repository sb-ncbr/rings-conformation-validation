#!/bin/bash

INPUT_DATA_FOLDER="$1"
OUTPUT_FOLDER="$2"

CYCLOPENTANE="cyclopentane"
CYCLOHEXANE="cyclohexane"
BENZENE="benzene"

python3 PrepareDataset.py -i "$INPUT_DATA_FOLDER" -o "$OUTPUT_FOLDER"
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "Error: PrepareDataset failed with exit code $exit_code"
    exit $exit_code
fi

for r in "$CYCLOPENTANE" "$CYCLOHEXANE"; do
    python3 FilterDataset.py -r "$r" -o "$OUTPUT_FOLDER" -i "$INPUT_DATA_FOLDER"
    exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "Error: FilterDataset failed with exit code $exit_code"
        exit $exit_code
    fi

    python3 CreateTemplates.py -r "$r" -o "$OUTPUT_FOLDER" -i "$INPUT_DATA_FOLDER"
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

python3 FilterDataset.py -r "$BENZENE" -o "$OUTPUT_FOLDER" -i "$INPUT_DATA_FOLDER"
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

CCP4="${INPUT_DATA_FOLDER}/ccp4"

python3 electron_density_coverage_analysis/main.py "$OUTPUT_FOLDER" "$CCP4"

python3 RingAnalysisResult.py -r "$CYCLOPENTANE" -o "$OUTPUT_FOLDER" -i "$INPUT_DATA_FOLDER"
python3 RingAnalysisResult.py -r "$CYCLOHEXANE" -o "$OUTPUT_FOLDER" -i "$INPUT_DATA_FOLDER"
python3 RingAnalysisResult.py -r "$BENZENE" -o "$OUTPUT_FOLDER" -i "$INPUT_DATA_FOLDER"