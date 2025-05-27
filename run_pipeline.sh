#!/bin/bash

# RNA3DB-2D-Structures End-to-End Pipeline
# This script orchestrates the entire workflow from downloading data to generating BPSEQ files.

# --- Default Configuration File ---
CONFIG_FILE="./config/pipeline_config.sh"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" # Gets the directory where the script is located

# --- Function to print messages ---
log_message() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] PIPELINE: $1" | tee -a "${PIPELINE_LOG_FILE_ABS}"
}

error_message() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] PIPELINE ERROR: $1" | tee -a "${PIPELINE_LOG_FILE_ABS}"
}

# --- Load Configuration ---
if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
    log_message "Loaded configuration from $CONFIG_FILE"
else
    error_message "Configuration file $CONFIG_FILE not found."
    error_message "Please copy 'config/pipeline_config_template.sh' to '$CONFIG_FILE' and edit it, or provide config via command line arguments."
    exit 1
fi

# --- Absolute path for log file ---
PIPELINE_LOG_FILE_ABS="${SCRIPT_DIR}/${PIPELINE_LOG_FILE}"
# Initialize log file
echo "Pipeline started on $(date)" > "${PIPELINE_LOG_FILE_ABS}"

# --- Create Output Directories ---
log_message "Creating output directories..."
mkdir -p "${SCRIPT_DIR}/data/${RAW_MMCIF_DIR}"
mkdir -p "${SCRIPT_DIR}/data/${ANALYSIS_OUTPUT_DIR}"
mkdir -p "${SCRIPT_DIR}/data/${BPSEQ_OUTPUT_DIR}"
log_message "Output directories created/ensured."

# --- Path to scripts ---
SCRIPTS_PATH="${SCRIPT_DIR}/scripts"

# --- Helper function to check command success ---
check_exit_status() {
    if [ $1 -ne 0 ]; then
        error_message "Error in $2 (Exit code: $1). See logs for details. Exiting pipeline."
        exit $1
    fi
}

# --- Start Pipeline ---
log_message "Starting RNA3DB-2D-Structures Pipeline..."
log_message "========================================="

# --- Step 1: Download Raw mmCIF files ---
log_message "[Step 1/6] Downloading raw mmCIF files..."
log_message "Input cluster.json: ${CLUSTER_JSON_PATH}"
log_message "Output directory for mmCIFs: ${SCRIPT_DIR}/${RAW_MMCIF_DIR}"
python "${SCRIPTS_PATH}/01_download_raw_mmcif.py" \
    "${CLUSTER_JSON_PATH}" \
    "${SCRIPT_DIR}/${RAW_MMCIF_DIR}" \
    --batch_script "${SCRIPT_DIR}/${BATCH_DOWNLOAD_SCRIPT}"
check_exit_status $? "Step 1: Download Raw mmCIF files"
log_message "Step 1 Complete."
log_message "-----------------------------------------"

# --- Step 2: Analyze RNA Chains ---
log_message "[Step 2/6] Analyzing RNA chains..."
log_message "Input mmCIF directory: ${SCRIPT_DIR}/${RAW_MMCIF_DIR}"
log_message "Output CSV: ${SCRIPT_DIR}/${RNA_CHAIN_ANALYSIS_CSV}"
log_message "Modifications cache: ${MODIFICATIONS_CACHE_PATH}"
log_message "Log file for analysis: ${SCRIPT_DIR}/${ANALYZE_LOG_FILE}"
python "${SCRIPTS_PATH}/02_analyze_rna_chains.py" \
    --input_dir "${SCRIPT_DIR}/${RAW_MMCIF_DIR}" \
    --output_csv "${SCRIPT_DIR}/${RNA_CHAIN_ANALYSIS_CSV}" \
    --mod_cache "${MODIFICATIONS_CACHE_PATH}" \ 
    --log_file "${SCRIPT_DIR}/${ANALYZE_LOG_FILE}" \
    ${ANALYZE_MAX_WORKERS:+--max_workers ${ANALYZE_MAX_WORKERS}}
check_exit_status $? "Step 2: Analyze RNA Chains"
log_message "Step 2 Complete."
log_message "-----------------------------------------"

# --- Step 3: Filter and Select Chains ---
log_message "[Step 3/6] Filtering and selecting representative chains..."
log_message "Input analysis CSV: ${SCRIPT_DIR}/${RNA_CHAIN_ANALYSIS_CSV}"
log_message "Input original cluster.json: ${CLUSTER_JSON_PATH}"
log_message "Output filtered/selected CSV: ${SCRIPT_DIR}/${FILTERED_SELECTED_CHAINS_CSV}"
log_message "Output selected_clusters.json for resplit: ${SCRIPT_DIR}/${SELECTED_CLUSTERS_JSON}"
python "${SCRIPTS_PATH}/03_filter_select_chains.py" \
    "${SCRIPT_DIR}/${RNA_CHAIN_ANALYSIS_CSV}" \
    "${CLUSTER_JSON_PATH}" \
    "${SCRIPT_DIR}/${FILTERED_SELECTED_CHAINS_CSV}" \
    "${SCRIPT_DIR}/${SELECTED_CLUSTERS_JSON}"
check_exit_status $? "Step 3: Filter and Select Chains"
log_message "Step 3 Complete."
log_message "-----------------------------------------"

# --- Step 4: Resplit Dataset ---
log_message "[Step 4/6] Resplitting dataset based on selected clusters..."
log_message "Input selected_clusters.json: ${SCRIPT_DIR}/${SELECTED_CLUSTERS_JSON}"
log_message "Output resplit_dataset.json: ${SCRIPT_DIR}/${RESPLIT_DATASET_JSON}"
bash "${SCRIPTS_PATH}/04_resplit_dataset.sh" \
    "${SCRIPT_DIR}/${SELECTED_CLUSTERS_JSON}" \
    "${SCRIPT_DIR}/${RESPLIT_DATASET_JSON}" \
    ${TRAIN_RATIO} \
    ${VALID_RATIO} \
    # Add logic for FORCE_ZERO_TEST if needed, e.g., by passing a flag
check_exit_status $? "Step 4: Resplit Dataset"
log_message "Step 4 Complete."
log_message "-----------------------------------------"

# --- Step 5: Assign Dataset Info ---
log_message "[Step 5/6] Assigning dataset information to selected chains..."
log_message "Input filtered/selected CSV: ${SCRIPT_DIR}/${FILTERED_SELECTED_CHAINS_CSV}"
log_message "Input resplit_dataset.json: ${SCRIPT_DIR}/${RESPLIT_DATASET_JSON}"
log_message "Output final all info CSV: ${SCRIPT_DIR}/${FINAL_ALL_INFO_CSV}"
log_message "Output final selected seqs CSV: ${SCRIPT_DIR}/${FINAL_SELECTED_SEQS_INFO_CSV}"
python "${SCRIPTS_PATH}/05_assign_dataset_info.py" \
    "${SCRIPT_DIR}/${FILTERED_SELECTED_CHAINS_CSV}" \
    "${SCRIPT_DIR}/${RESPLIT_DATASET_JSON}" \
    "${SCRIPT_DIR}/${FINAL_ALL_INFO_CSV}" \
    "${SCRIPT_DIR}/${FINAL_SELECTED_SEQS_INFO_CSV}"
check_exit_status $? "Step 5: Assign Dataset Info"
log_message "Step 5 Complete."
log_message "-----------------------------------------"

# --- Step 6: Generate BPSEQ files ---
log_message "[Step 6/6] Generating BPSEQ files..."
log_message "Input final selected seqs CSV: ${SCRIPT_DIR}/${FINAL_SELECTED_SEQS_INFO_CSV}"
log_message "Output BPSEQ directory: ${SCRIPT_DIR}/${BPSEQ_OUTPUT_DIR}"
GENERATE_LESS_600NT_FLAG=""
if [ "$GENERATE_LESS_600NT_BPSEQ" = true ]; then
    GENERATE_LESS_600NT_FLAG="--generate_less_600nt"
fi
python "${SCRIPTS_PATH}/06_generate_bpseq.py" \
    "${SCRIPT_DIR}/${FINAL_SELECTED_SEQS_INFO_CSV}" \
    "${SCRIPT_DIR}/${BPSEQ_OUTPUT_DIR}" \
    ${GENERATE_LESS_600NT_FLAG} \
    --less_600nt_subdir_name "${BPSEQ_LESS_600NT_SUBDIR_NAME}"
check_exit_status $? "Step 6: Generate BPSEQ files"
log_message "Step 6 Complete."
log_message "-----------------------------------------"

log_message "========================================="
log_message "Pipeline finished successfully!"
log_message "Final BPSEQ files are in: ${SCRIPT_DIR}/${BPSEQ_OUTPUT_DIR}"
log_message "Final selected sequences information CSV: ${SCRIPT_DIR}/${FINAL_SELECTED_SEQS_INFO_CSV}"
log_message "Full pipeline log: ${PIPELINE_LOG_FILE_ABS}"

exit 0