#!/bin/bash

# This script re-splits an RNA3DB cluster JSON file (typically one
# containing only selected representative clusters) into train, validation,
# and test sets using the rna3db library's split functionality.
#
# 此脚本使用 rna3db 库的 split 功能，将一个 RNA3DB 集群 JSON 文件
# （通常是仅包含选定的代表性集群的文件）重新划分为训练集、验证集和测试集。

# --- Script Safety ---
set -e # Exit immediately if a command exits with a non-zero status.
set -u # Treat unset variables as an error when substituting.
set -o pipefail # Return value of a pipeline is the value of the last command to exit with a non-zero status

# --- Argument Parsing ---
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <input_selected_clusters.json> <output_resplit.json> <train_ratio> <valid_ratio> [--force_zero_test]"
    echo "  <input_selected_clusters.json>: Path to the JSON file containing clusters to be split."
    echo "                                    (e.g., output from 03_filter_select_chains.py)"
    echo "  <output_resplit.json>: Path for the output JSON file with dataset splits."
    echo "  <train_ratio>: Proportion for the training set (e.g., 0.6)."
    echo "  <valid_ratio>: Proportion for the validation set (e.g., 0.2)."
    echo "  [--force_zero_test]: Optional. If provided, forces component_0 (if present) into the test set."
    exit 1
fi

INPUT_JSON="$1"
OUTPUT_JSON="$2"
TRAIN_RATIO="$3"
VALID_RATIO="$4"
FORCE_ZERO_TEST_FLAG=""

# Check for optional flag
if [[ "${5:-}" == "--force_zero_test" ]]; then
    FORCE_ZERO_TEST_FLAG="--force_zero_test"
    echo "INFO: --force_zero_test flag is set. Component_0 will be forced into the test set if present."
fi

# --- Validation ---
if [ ! -f "$INPUT_JSON" ]; then
    echo "ERROR: Input JSON file not found: $INPUT_JSON"
    exit 1
fi

# Validate ratios (basic check, rna3db will do more thorough checks)
if (( $(echo "$TRAIN_RATIO < 0 || $TRAIN_RATIO > 1" | bc -l) )); then
    echo "ERROR: Train ratio must be between 0 and 1."
    exit 1
fi
if (( $(echo "$VALID_RATIO < 0 || $VALID_RATIO > 1" | bc -l) )); then
    echo "ERROR: Validation ratio must be between 0 and 1."
    exit 1
fi
TOTAL_RATIO=$(echo "$TRAIN_RATIO + $VALID_RATIO" | bc -l)
if (( $(echo "$TOTAL_RATIO > 1" | bc -l) )); then
    echo "ERROR: Sum of train_ratio and valid_ratio cannot exceed 1."
    exit 1
fi

# --- Ensure output directory exists ---
OUTPUT_DIR=$(dirname "$OUTPUT_JSON")
mkdir -p "$OUTPUT_DIR"
if [ $? -ne 0 ]; then
    echo "ERROR: Could not create output directory: $OUTPUT_DIR"
    exit 1
fi

# --- Run rna3db split command ---
# The rna3db library needs to be installed and accessible in the environment.
# The command structure is: python -m rna3db split <input> <output> --train_ratio X --valid_ratio Y [--force_zero_test]
echo "INFO: Running rna3db split..."
echo "  Input: $INPUT_JSON"
echo "  Output: $OUTPUT_JSON"
echo "  Train Ratio: $TRAIN_RATIO"
echo "  Valid Ratio: $VALID_RATIO"
if [ -n "$FORCE_ZERO_TEST_FLAG" ]; then
    echo "  Force Zero Test: Enabled"
fi

# Construct the command
COMMAND=(python -m rna3db split \
    "$INPUT_JSON" \
    "$OUTPUT_JSON" \
    --train_ratio "$TRAIN_RATIO" \
    --valid_ratio "$VALID_RATIO")

if [ -n "$FORCE_ZERO_TEST_FLAG" ]; then
    COMMAND+=("$FORCE_ZERO_TEST_FLAG")
fi

# Execute the command
"${COMMAND[@]}"

# Check execution result
if [ $? -ne 0 ]; then
    echo "ERROR: rna3db split command failed. Check rna3db logs or output for details."
    exit 1
else
    echo "INFO: Successfully generated resplit dataset: $OUTPUT_JSON"
fi

exit 0