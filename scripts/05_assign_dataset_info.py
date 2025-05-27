#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to assign dataset partition information (train/valid/test) from a
resplit JSON file to a DataFrame of filtered and selected RNA chains.

This script performs operations similar to step 9 in the original
RDL1.6.7-1.5_valid_chain_selection.ipynb notebook.
"""

import pandas as pd
import json
import argparse
import sys
from pathlib import Path
import logging

# --- Global Logger ---
logger = logging.getLogger(__name__)

# --- Logging Setup ---
def setup_logging(log_level_str: str = "INFO", log_filepath: Optional[Path] = None):
    log_level = getattr(logging, log_level_str.upper(), logging.INFO)
    log_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s')
    
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    if root_logger.hasHandlers():
        root_logger.handlers.clear()

    if log_filepath:
        try:
            file_handler = logging.FileHandler(log_filepath, mode='w', encoding='utf-8')
            file_handler.setFormatter(log_formatter)
            file_handler.setLevel(log_level)
            root_logger.addHandler(file_handler)
        except Exception as e:
            print(f"Error: Could not set up log file handler at {log_filepath}: {e}")

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_formatter)
    console_handler.setLevel(log_level)
    root_logger.addHandler(console_handler)
    logger.info(f"Logging setup complete. Level: {log_level_str}{', File: ' + str(log_filepath) if log_filepath else ''}")


def parse_resplit_json_to_dataset_map(resplit_json_path: Path) -> dict:
    """
    Parses the resplit JSON (output of 04_resplit_dataset.sh) and creates a
    mapping from 'repr_pdb_chain' to its assigned dataset (train_set, valid_set, test_set).

    解析 resplit JSON（04_resplit_dataset.sh 的输出）并创建一个从 'repr_pdb_chain'
    到其分配的数据集（train_set, valid_set, test_set）的映射。

    Args:
        resplit_json_path (Path): Path to the resplit JSON file.
                                  resplit JSON 文件的路径。
    Returns:
        dict: A dictionary mapping repr_pdb_chain (str) to dataset_name (str).
              一个将 repr_pdb_chain (str) 映射到 dataset_name (str) 的字典。
    """
    logger.info(f"Parsing resplit dataset JSON from: {resplit_json_path}")
    with open(resplit_json_path, 'r') as f:
        resplit_data = json.load(f)
    
    repr_to_dataset_map = {}
    # Structure of resplit_data:
    # {
    #   "train_set": {"component_X": {"repr_pdb_A": {...}, "repr_pdb_B": {...}}, ...},
    #   "valid_set": {"component_Y": {"repr_pdb_C": {...}}, ...},
    #   "test_set":  {"component_Z": {"repr_pdb_D": {...}}, ...}
    # }
    for dataset_name, components in resplit_data.items():
        for component_id, repr_pdb_clusters in components.items():
            for repr_pdb_chain_id in repr_pdb_clusters.keys():
                if repr_pdb_chain_id in repr_to_dataset_map:
                    logger.warning(f"Representative PDB chain '{repr_pdb_chain_id}' found in multiple datasets ('{repr_to_dataset_map[repr_pdb_chain_id]}' and '{dataset_name}'). This should not happen with rna3db split. Using the first encountered: '{repr_to_dataset_map[repr_pdb_chain_id]}'.")
                else:
                    repr_to_dataset_map[repr_pdb_chain_id] = dataset_name
    
    logger.info(f"Created map for {len(repr_to_dataset_map)} representative PDB chains to their datasets.")
    return repr_to_dataset_map


def main():
    parser = argparse.ArgumentParser(description="Assign dataset partition info to filtered/selected RNA chains.")
    parser.add_argument("filtered_selected_csv", type=Path, help="Path to the CSV from 03_filter_select_chains.py (output_filtered_selected_csv).")
    parser.add_argument("resplit_json", type=Path, help="Path to the resplit JSON file from 04_resplit_dataset.sh.")
    parser.add_argument("output_final_all_info_csv", type=Path, help="Path to save the final CSV containing ALL chains with dataset assignments.")
    parser.add_argument("output_final_selected_seqs_csv", type=Path, help="Path to save the final CSV containing ONLY SELECTED chains with dataset assignments (for BPSEQ generation).")
    parser.add_argument("--log_file", type=Path, default=None, help="Optional path to a log file.")
    parser.add_argument("--log_level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="Logging level.")

    args = parser.parse_args()

    setup_logging(log_level_str=args.log_level, log_filepath=args.log_file)

    # --- 1. Load filtered and selected chains CSV ---
    logger.info(f"Loading filtered and selected chains data from: {args.filtered_selected_csv}")
    try:
        df_chains = pd.read_csv(args.filtered_selected_csv)
    except FileNotFoundError:
        logger.error(f"Filtered/selected chains CSV file not found: {args.filtered_selected_csv}")
        sys.exit(1)
    logger.info(f"Loaded {len(df_chains)} chains.")
    
    if 'repr_pdb_chain' not in df_chains.columns:
        logger.error(f"Critical column 'repr_pdb_chain' not found in {args.filtered_selected_csv}. This column is needed to map dataset assignments.")
        sys.exit(1)
    if 'is_selected' not in df_chains.columns:
        logger.warning(f"Column 'is_selected' not found in {args.filtered_selected_csv}. Will assume all chains are to be processed if filtering for selected ones later.")
        # Add a default 'is_selected' = 1 if it's missing, so the script doesn't break,
        # but this indicates a potential issue in the upstream script (03).
        df_chains['is_selected'] = 1


    # --- 2. Parse resplit JSON to get dataset assignments for representative PDBs ---
    repr_to_dataset_map = parse_resplit_json_to_dataset_map(args.resplit_json)
    if not repr_to_dataset_map:
        logger.warning(f"The resplit JSON ({args.resplit_json}) resulted in an empty map. No dataset assignments can be made.")
        # Create an empty 'dataset' column and save.
        df_chains['dataset'] = "unassigned"
    else:
        # --- 3. Map dataset assignment to each chain in df_chains ---
        # Each chain in df_chains belongs to a 'repr_pdb_chain'.
        # We use this 'repr_pdb_chain' to look up its dataset from the map.
        logger.info("Mapping dataset assignments to chains...")
        df_chains['dataset'] = df_chains['repr_pdb_chain'].map(repr_to_dataset_map)
        
        unassigned_count = df_chains['dataset'].isna().sum()
        if unassigned_count > 0:
            logger.warning(f"{unassigned_count} chains could not be assigned to a dataset. Their 'repr_pdb_chain' might not be in the resplit JSON. They will be marked as 'unassigned'.")
            df_chains['dataset'].fillna('unassigned', inplace=True)

    # --- 4. Save the DataFrame with all chains and their dataset assignments ---
    df_chains.to_csv(args.output_final_all_info_csv, index=False)
    logger.info(f"Final CSV with all chains and dataset assignments saved to: {args.output_final_all_info_csv}")

    # --- 5. Filter for only selected chains and save ---
    logger.info("Filtering for 'is_selected == 1' chains for BPSEQ generation input...")
    df_final_selected_seqs = df_chains[df_chains['is_selected'] == 1].copy()
    
    # Sanity check: ensure all selected sequences were assigned a dataset (unless map was empty)
    if repr_to_dataset_map: # Only check if assignments were expected
        unassigned_selected_count = df_final_selected_seqs[df_final_selected_seqs['dataset'] == 'unassigned'].shape[0]
        if unassigned_selected_count > 0:
             logger.error(f"CRITICAL: {unassigned_selected_count} chains marked as 'is_selected == 1' were NOT assigned a dataset. This indicates an inconsistency between selected representatives and the resplit JSON content.")
             # This could happen if a 'repr_pdb_chain' that had a selected member was somehow excluded from the input to 04_resplit_dataset.sh
        else:
            logger.info(f"All {len(df_final_selected_seqs)} selected chains were successfully assigned a dataset (or marked 'unassigned' if resplit map was empty).")


    df_final_selected_seqs.to_csv(args.output_final_selected_seqs_csv, index=False)
    logger.info(f"Final CSV with ONLY SELECTED chains and dataset assignments saved to: {args.output_final_selected_seqs_csv}")
    logger.info(f"This file contains {len(df_final_selected_seqs)} sequences and will be used for BPSEQ generation.")
    
    logger.info("Script finished successfully.")

if __name__ == "__main__":
    main()