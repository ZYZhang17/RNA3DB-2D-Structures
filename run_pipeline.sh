#!/bin/bash

# -----------------------------------------------------------------------------
# RNA3DB-2D-Structures Pipeline Runner
#
# This script executes the complete pipeline to:
# 1. Download mmCIF files for PDB IDs in an RNA3DB cluster.json
# 2. Analyze RNA chains and extract secondary structures
# 3. Filter and select representative chains
# 4. Re-split the dataset into train/valid/test sets
# 5. Assign dataset information to selected chains
# 6. Generate BPSEQ files for the selected chains
# 7. Convert BPSEQ files to pickle format for model training
# -----------------------------------------------------------------------------

# === Script Safety ===
set -e  # Exit immediately if a command exits with non-zero status
set -u  # Treat unset variables as an error
set -o pipefail  # Return value of a pipeline is the value of the last command to exit with a non-zero status

# === Load Configuration ===
# Check if pipeline_config.sh exists
if [ ! -f "./pipeline_config.sh" ]; then
    echo "Error: pipeline_config.sh not found!"
    echo "Please create it by copying and editing pipeline_config_template.sh:"
    echo "  cp pipeline_config_template.sh pipeline_config.sh"
    exit 1
fi

# Source the configuration file
source ./pipeline_config.sh

# === Setup Logging ===
# Create analysis output directory if it doesn't exist
mkdir -p "${ANALYSIS_OUTPUT_DIR}"

# Start the pipeline log
echo "==== RNA3DB-2D-Structures Pipeline Started: $(date) ====" > "${PIPELINE_LOG_FILE}"
echo "Using configuration from: ./pipeline_config.sh" >> "${PIPELINE_LOG_FILE}"

# === Function to log commands ===
log_cmd() {
    echo -e "\n>>> Running: $*\n" | tee -a "${PIPELINE_LOG_FILE}"
    "$@" 2>&1 | tee -a "${PIPELINE_LOG_FILE}"
    return ${PIPESTATUS[0]}
}

# === Step 1: Download Raw mmCIF Files ===
echo "=== Step 1: Downloading Raw mmCIF Files ===" | tee -a "${PIPELINE_LOG_FILE}"
mkdir -p "${RAW_MMCIF_DIR}"

log_cmd python3 ./scripts/01_download_raw_mmcif.py \
    "${CLUSTER_JSON_PATH}" \
    "${RAW_MMCIF_DIR}" \
    --batch_script "${BATCH_DOWNLOAD_SCRIPT}"

# === Step 2: Analyze RNA Chains ===
echo "=== Step 2: Analyzing RNA Chains ===" | tee -a "${PIPELINE_LOG_FILE}"

log_cmd python3 ./scripts/02_analyze_rna_chains.py \
    --input_dir "${RAW_MMCIF_DIR}" \
    --output_csv "${RNA_CHAIN_ANALYSIS_CSV}" \
    --mod_cache "${MODIFICATIONS_CACHE_PATH}" \
    --log_file "${ANALYZE_LOG_FILE}" \
    --max_workers "${ANALYZE_MAX_WORKERS}"

# === Step 3: Filter and Select Chains ===
echo "=== Step 3: Filtering and Selecting Chains ===" | tee -a "${PIPELINE_LOG_FILE}"

log_cmd python3 ./scripts/03_filter_select_chains.py \
    --analysis_csv "${RNA_CHAIN_ANALYSIS_CSV}" \
    --cluster_json "${CLUSTER_JSON_PATH}" \
    --output_filtered_selected_csv "${FILTERED_SELECTED_CHAINS_CSV}" \
    --output_selected_clusters_json "${SELECTED_CLUSTERS_JSON}" \
    --min_valid_proportion "${FILTER_MIN_VALID_PROPORTION}" \
    --random_seed "${FILTER_RANDOM_SEED}" \
    --log_file "${FILTER_SELECT_LOG_FILE}" \
    --log_level "${DEFAULT_LOG_LEVEL}"

# === Step 4: Resplit Dataset ===
echo "=== Step 4: Resplitting Dataset ===" | tee -a "${PIPELINE_LOG_FILE}"

# Prepare force_zero_test flag if needed
FORCE_ZERO_TEST_FLAG=""
if [ "${FORCE_ZERO_TEST}" = "true" ]; then
    FORCE_ZERO_TEST_FLAG="--force_zero_test"
fi

log_cmd bash ./scripts/04_resplit_dataset.sh \
    "${SELECTED_CLUSTERS_JSON}" \
    "${RESPLIT_DATASET_JSON}" \
    "${TRAIN_RATIO}" \
    "${VALID_RATIO}" \
    ${FORCE_ZERO_TEST_FLAG}

# === Step 5: Assign Dataset Info ===
echo "=== Step 5: Assigning Dataset Information ===" | tee -a "${PIPELINE_LOG_FILE}"

log_cmd python3 ./scripts/05_assign_dataset_info.py \
    --filtered_selected_csv "${FILTERED_SELECTED_CHAINS_CSV}" \
    --resplit_json "${RESPLIT_DATASET_JSON}" \
    --output_final_all_info_csv "${FINAL_ALL_INFO_CSV}" \
    --output_final_selected_dataset_csv "${FINAL_SELECTED_SEQS_INFO_CSV}" \
    --log_file "${ASSIGN_DATASET_LOG_FILE}" \
    --log_level "${DEFAULT_LOG_LEVEL}"

# === Step 6: Generate BPSEQ Files ===
echo "=== Step 6: Generating BPSEQ Files ===" | tee -a "${PIPELINE_LOG_FILE}"

# Prepare less_600nt flag if needed
GENERATE_LESS_600NT_FLAG=""
if [ "${GENERATE_LESS_600NT_BPSEQ}" = "true" ]; then
    GENERATE_LESS_600NT_FLAG="--generate_less_600nt"
fi

log_cmd python3 ./scripts/06_generate_bpseq.py \
    --input_selected_seqs_csv "${FINAL_SELECTED_SEQS_INFO_CSV}" \
    --output_bpseq_dir "${BPSEQ_OUTPUT_DIR}" \
    ${GENERATE_LESS_600NT_FLAG} \
    --less_600nt_subdir_name "${BPSEQ_LESS_600NT_SUBDIR_NAME}" \
    --log_file "${GENERATE_BPSEQ_LOG_FILE}" \
    --log_level "${DEFAULT_LOG_LEVEL}"

# === Step 7: Convert BPSEQ to Pickle ===
echo "=== Step 7: Converting BPSEQ to Pickle Format ===" | tee -a "${PIPELINE_LOG_FILE}"

# Create the pickle output directory
mkdir -p "${PICKLE_OUTPUT_DIR}"

# Process each dataset
for dataset in ${PICKLE_DATASETS}; do
    echo "Processing dataset: ${dataset}" | tee -a "${PIPELINE_LOG_FILE}"
    
    # Construct the BPSEQ directory for this dataset
    bpseq_dataset_dir="${BPSEQ_INPUT_BASE_DIR}/${dataset}"
    
    # Construct the output pickle file path
    output_pickle="${PICKLE_OUTPUT_DIR}/${dataset}.pickle"
    
    # Check if BPSEQ directory for this dataset exists
    if [ ! -d "${bpseq_dataset_dir}" ]; then
        echo "Warning: BPSEQ directory for ${dataset} not found: ${bpseq_dataset_dir}" | tee -a "${PIPELINE_LOG_FILE}"
        continue
    fi
    
    log_cmd python3 ./scripts/07_bpseq_to_pickle.py \
        --bpseq_dir "${bpseq_dataset_dir}" \
        --output_pickle "${output_pickle}" \
        --split_json "${PICKLE_SPLIT_JSON_PATH}" \
        --dataset_name "${dataset}" \
        --random_seed "${PICKLE_RANDOM_SEED}" \
        --max_length "${PICKLE_MAX_LENGTH}" \
        --log_file "${BPSEQ_TO_PICKLE_LOG_FILE}" \
        --log_level "${DEFAULT_LOG_LEVEL}"
done

# === Pipeline Complete ===
echo "==== RNA3DB-2D-Structures Pipeline Completed: $(date) ====" | tee -a "${PIPELINE_LOG_FILE}"
echo "All steps executed successfully."
echo "Pipeline log: ${PIPELINE_LOG_FILE}"