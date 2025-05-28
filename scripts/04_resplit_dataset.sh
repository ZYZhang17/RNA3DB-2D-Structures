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

# --- Default Values ---
DEFAULT_INPUT_JSON="./data/analysis_output/03_selected_clusters.json"
DEFAULT_OUTPUT_JSON="./data/analysis_output/04_resplit.json"
DEFAULT_TRAIN_RATIO="0.6"
DEFAULT_VALID_RATIO="0.2"
DEFAULT_FORCE_ZERO_TEST="True"

# --- Argument Parsing ---
# 使用提供的参数或默认值
INPUT_JSON="${1:-$DEFAULT_INPUT_JSON}"
OUTPUT_JSON="${2:-$DEFAULT_OUTPUT_JSON}"
TRAIN_RATIO="${3:-$DEFAULT_TRAIN_RATIO}"
VALID_RATIO="${4:-$DEFAULT_VALID_RATIO}"

# 第5个参数如果是--force_zero_test或未提供且默认为True，则设置标志
if [[ "${5:-}" == "--force_zero_test" || ($# -lt 5 && "$DEFAULT_FORCE_ZERO_TEST" == "True") ]]; then
    FORCE_ZERO_TEST_FLAG="--force_zero_test"
    echo "INFO: --force_zero_test flag is set. Component_0 will be forced into the test set if present."
else
    FORCE_ZERO_TEST_FLAG=""
fi

# 显示使用说明（如果使用-h或--help）
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    echo "Usage: $0 [<input_json>] [<output_json>] [<train_ratio>] [<valid_ratio>] [--force_zero_test]"
    echo "  <input_json>: Path to the JSON file containing clusters to be split."
    echo "                 Default: $DEFAULT_INPUT_JSON"
    echo "  <output_json>: Path for the output JSON file with dataset splits."
    echo "                 Default: $DEFAULT_OUTPUT_JSON"
    echo "  <train_ratio>: Proportion for the training set."
    echo "                 Default: $DEFAULT_TRAIN_RATIO"
    echo "  <valid_ratio>: Proportion for the validation set."
    echo "                 Default: $DEFAULT_VALID_RATIO"
    echo "  --force_zero_test: If provided, forces component_0 (if present) into the test set."
    echo "                     Default: $(if [[ "$DEFAULT_FORCE_ZERO_TEST" == "True" ]]; then echo "Enabled"; else echo "Disabled"; fi)"
    exit 0
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

# --- Run local split command ---
echo "INFO: Running local split function..."
echo "  Input: $INPUT_JSON"
echo "  Output: $OUTPUT_JSON"
echo "  Train Ratio: $TRAIN_RATIO"
echo "  Valid Ratio: $VALID_RATIO"
if [ -n "$FORCE_ZERO_TEST_FLAG" ]; then
    echo "  Force Zero Test: Enabled"
fi

# Calculate test ratio
TEST_RATIO=$(python -c "print(1.0 - $TRAIN_RATIO - $VALID_RATIO)")

# Create Python script to call split function
PYTHON_SCRIPT="
import sys
from pathlib import Path

# Add scripts directory to Python path
scripts_dir = Path('$0').parent.absolute()
sys.path.insert(0, str(scripts_dir))

from utils import split

# Call split function
split(
    input_path='$INPUT_JSON',
    output_path='$OUTPUT_JSON',
    splits=[$TRAIN_RATIO, $VALID_RATIO, $TEST_RATIO],
    split_names=['train_set', 'valid_set', 'test_set'],
    shuffle=False,
    force_zero_last=$(if [ -n "$FORCE_ZERO_TEST_FLAG" ]; then echo "True"; else echo "False"; fi)
)
"

# Execute the Python script
python -c "$PYTHON_SCRIPT"

# Check execution result
if [ $? -ne 0 ]; then
    echo "ERROR: split command failed. Check logs for details."
    exit 1
else
    echo "INFO: Successfully generated resplit dataset: $OUTPUT_JSON"
fi

exit 0