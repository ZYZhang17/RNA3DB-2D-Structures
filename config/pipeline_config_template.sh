#!/bin/bash

# -----------------------------------------------------------------------------
# Configuration for RNA3DB-2D-Structures Pipeline
#
# Instructions:
# 1. Copy this file to 'pipeline_config.sh' in the same directory:
#    cp pipeline_config_template.sh pipeline_config.sh
# 2. Edit 'pipeline_config.sh' with your specific paths and settings.
#    This file (pipeline_config.sh) will be ignored by git.
# -----------------------------------------------------------------------------

# === Input Data ===
# Path to your cluster.json file from an RNA3DB release
# (e.g., from rna3db-jsons.tar.gz)
CLUSTER_JSON_PATH="./data/RNA3DB_2024-12-04-release_cluster.json"

# Path to the modifications_cache.json file
# The project provides one in ./data/modifications_cache.json
MODIFICATIONS_CACHE_PATH="./utils/modifications_cache.json"


# === External Tools ===
# Path to the batch_download.sh script
# The project provides one in ./utils/batch_download.sh
BATCH_DOWNLOAD_SCRIPT="./utils/batch_download.sh"


# === Output Directories and Files ===
# Directory where raw mmCIF files will be downloaded
RAW_MMCIF_DIR="./data/downloaded_mmcif"

# Directory for analysis outputs (CSVs, intermediate JSONs)
ANALYSIS_OUTPUT_DIR="./data/analysis_output"

# Directory for final BPSEQ files
BPSEQ_OUTPUT_DIR="./bpseq_files"

# --- Specific output file names (can be kept as default) ---
# CSV from RNA chain analysis (Step 2 output)
RNA_CHAIN_ANALYSIS_CSV="${ANALYSIS_OUTPUT_DIR}/02_rna_chain_analysis.csv"

# JSON file containing selected clusters for resplitting (Step 3 output)
SELECTED_CLUSTERS_JSON="${ANALYSIS_OUTPUT_DIR}/selected_clusters.json"

# CSV with all filtered chains and selection marks (Step 3 output)
FILTERED_SELECTED_CHAINS_CSV="${ANALYSIS_OUTPUT_DIR}/filtered_selected_chains.csv"

# JSON file from dataset resplitting (Step 4 output)
RESPLIT_DATASET_JSON="${ANALYSIS_OUTPUT_DIR}/resplit_dataset.json"

# Final CSV with all information including dataset assignment (Step 5 output)
FINAL_ALL_INFO_CSV="${ANALYSIS_OUTPUT_DIR}/final_all_info.csv"

# Final CSV with only selected sequences and dataset assignment (Step 5 output, input for BPSEQ)
FINAL_SELECTED_SEQS_INFO_CSV="${ANALYSIS_OUTPUT_DIR}/final_selected_seqs_info.csv"

# BPSEQ subdirectory for sequences < 600nt (optional, handled in 06_generate_bpseq.py)
BPSEQ_LESS_600NT_SUBDIR_NAME="bpseq_less600nt" # This will be inside BPSEQ_OUTPUT_DIR


# === Script Specific Settings ===
# For 02_analyze_rna_chains.py
# Number of parallel workers for mmCIF analysis (set to None for os.cpu_count())
ANALYZE_MAX_WORKERS=None # Example: 8 or None
ANALYZE_LOG_FILE="${ANALYSIS_OUTPUT_DIR}/02_analyze_rna_chains.log"

# For 04_resplit_dataset.sh (RNA3DB split ratios)
TRAIN_RATIO=0.6
VALID_RATIO=0.2
# TEST_RATIO is implicitly 1 - TRAIN_RATIO - VALID_RATIO
# FORCE_ZERO_TEST can be true or false (passed as a flag to the script)

# For 06_generate_bpseq.py
# Whether to also generate a separate set of BPSEQ files for sequences < 600nt
GENERATE_LESS_600NT_BPSEQ=true


# === Logging ===
PIPELINE_LOG_FILE="./pipeline.log"

# End of Configuration