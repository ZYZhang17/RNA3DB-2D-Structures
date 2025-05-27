#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to filter RNA chains from analysis results, select the best representative
for each RNA3DB cluster, and prepare a new cluster file for dataset splitting.

This script performs operations similar to steps 1-7 in the original
RDL1.6.7-1.5_valid_chain_selection.ipynb notebook.
"""

import pandas as pd
import json
import numpy as np
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


def calculate_seq_valid_proportion(df_analysis_results: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the proportion of the RNApolis raw sequence length to the canonical sequence length.
    Handles potential division by zero or NaN lengths by returning 0 for such cases.

    计算 RNApolis 原始序列长度与权威序列长度的比例。
    通过在此类情况下返回0来处理潜在的除零或NaN长度。

    Args:
        df_analysis_results (pd.DataFrame): DataFrame from 02_analyze_rna_chains.py output.
                                            来自 02_analyze_rna_chains.py 输出的 DataFrame。
    Returns:
        pd.DataFrame: DataFrame with a new 'seq_valid_proportion' column.
                      带有新的 'seq_valid_proportion' 列的 DataFrame。
    """
    logger.info("Calculating sequence valid proportion...")
    
    # Ensure lengths are numeric, fill NaNs with 0 to avoid errors in len() or division
    # Sequence_RNApolis_Raw might be "N/A" if RNApolis failed for a chain
    rnapolis_lengths = df_analysis_results['Sequence_RNApolis_Raw'].apply(lambda x: len(x) if isinstance(x, str) and x != "N/A" else 0)
    canonical_lengths = df_analysis_results['Sequence_Canonical'].apply(lambda x: len(x) if isinstance(x, str) and x != "N/A" else 0)

    # Calculate proportion, handling division by zero
    df_analysis_results['seq_valid_proportion'] = np.where(
        canonical_lengths > 0, 
        rnapolis_lengths / canonical_lengths, 
        0 # Set to 0 if canonical length is 0 or N/A
    )
    logger.info("Sequence valid proportion calculation complete.")
    return df_analysis_results


def rna3db_cluster_json_to_df(cluster_json_path: Path) -> pd.DataFrame:
    """
    Reads an RNA3DB cluster.json file and flattens it into a Pandas DataFrame.
    读取 RNA3DB cluster.json 文件并将其扁平化为 Pandas DataFrame。

    Args:
        cluster_json_path (Path): Path to the cluster.json file.
                                  cluster.json 文件路径。
    Returns:
        pd.DataFrame: DataFrame with columns: component, repr_pdb_chain, pdb_chain_id,
                      and other metadata from the JSON.
                      包含列：component, repr_pdb_chain, pdb_chain_id 以及来自 JSON 的其他元数据。
    """
    logger.info(f"Loading and parsing RNA3DB cluster.json from: {cluster_json_path}")
    with open(cluster_json_path, 'r') as f:
        data = json.load(f)
    
    rows = []
    # cluster.json structure: {component_id: {repr_pdb_id: {actual_pdb_chain_id: metadata}}}
    for component_id, component_data in data.items():
        for repr_pdb_id, repr_cluster_data in component_data.items():
            for actual_pdb_chain_id, metas in repr_cluster_data.items():
                row = {
                    "component": component_id,
                    "repr_pdb_chain": repr_pdb_id, # This is the representative of the cluster
                    "pdb_chain_id_from_cluster": actual_pdb_chain_id # This is the specific chain
                }
                row.update(metas) # Add all metadata like release_date, resolution, etc.
                rows.append(row)
    
    json_df = pd.DataFrame(rows)
    logger.info(f"Parsed cluster.json into DataFrame with {len(json_df)} rows.")
    return json_df


def select_representative_chain(df_group: pd.DataFrame, random_seed: int = 66) -> pd.Series:
    """
    Selects the best representative chain from a group based on criteria:
    1. Lowest (best) resolution.
    2. Highest 'seq_valid_proportion'.
    3. Random choice if ties persist.

    从一组中选择最佳代表链，基于以下标准：
    1. 最低（最佳）分辨率。
    2. 最高的 'seq_valid_proportion'。
    3. 如果仍然存在平局，则随机选择。

    Args:
        df_group (pd.DataFrame): DataFrame group (all chains belonging to one 'repr_pdb_chain').
                                 DataFrame 组（属于一个 'repr_pdb_chain' 的所有链）。
        random_seed (int): Seed for random selection.
                           随机选择的种子。
    Returns:
        pd.Series: The row of the selected representative chain.
                   所选代表链的行。
    """
    # Ensure resolution is numeric, coercing errors to NaN, then fill NaN with a high value for min()
    df_group['resolution_numeric'] = pd.to_numeric(df_group['resolution'], errors='coerce').fillna(float('inf'))

    # 1. Resolution: Lower is better
    min_res = df_group['resolution_numeric'].min()
    df1 = df_group[df_group['resolution_numeric'] == min_res]
    
    # 2. seq_valid_proportion: Higher is better
    max_seq_prop = df1['seq_valid_proportion'].max()
    df2 = df1[df1['seq_valid_proportion'] == max_seq_prop]
    
    # 3. If still multiple, random choice
    if len(df2) > 1:
        # Use a temporary random state for this choice to not affect global numpy random state
        # if this function is called multiple times (e.g. in apply)
        temp_rng = np.random.RandomState(random_seed)
        selected_row = df2.sample(n=1, random_state=temp_rng).iloc[0]
    elif len(df2) == 1:
        selected_row = df2.iloc[0]
    else: # Should not happen if df_group is not empty
        logger.warning(f"No chain selected for group. Group was: {df_group['repr_pdb_chain'].iloc[0] if not df_group.empty else 'Empty Group'}")
        return pd.Series(dtype='object') # Return empty series
        
    return selected_row


def filter_rna3db_json_by_repr_ids(input_json_path: Path, selected_repr_pdb_ids: set, output_json_path: Path):
    """
    Filters the original cluster.json to keep only the specified representative PDB chain clusters.
    过滤原始 cluster.json，仅保留指定的代表性 PDB 链簇。

    Args:
        input_json_path (Path): Path to the original cluster.json.
                                原始 cluster.json 的路径。
        selected_repr_pdb_ids (set): Set of 'repr_pdb_chain' IDs to keep.
                                     要保留的 'repr_pdb_chain' ID 集合。
        output_json_path (Path): Path to save the filtered JSON.
                                 保存过滤后的 JSON 的路径。
    """
    logger.info(f"Filtering original cluster.json based on {len(selected_repr_pdb_ids)} selected representative PDB chain IDs.")
    with open(input_json_path, 'r') as f:
        data = json.load(f)
    
    filtered_data = {}
    kept_clusters_count = 0
    for component_id, component_data in data.items():
        filtered_pdbs_in_component = {}
        for repr_pdb_id, repr_cluster_data in component_data.items():
            if repr_pdb_id in selected_repr_pdb_ids:
                filtered_pdbs_in_component[repr_pdb_id] = repr_cluster_data
                kept_clusters_count +=1
        
        if filtered_pdbs_in_component:
            filtered_data[component_id] = filtered_pdbs_in_component
            
    with open(output_json_path, 'w') as f:
        json.dump(filtered_data, f, indent=2)
    logger.info(f"Filtered cluster.json saved to: {output_json_path}. Kept {kept_clusters_count} representative clusters.")


def main():
    parser = argparse.ArgumentParser(description="Filter RNA chains, select representatives, and prepare for dataset splitting.")
    parser.add_argument("analysis_csv", type=Path, help="Path to the CSV file from 02_analyze_rna_chains.py.")
    parser.add_argument("cluster_json", type=Path, help="Path to the original RNA3DB cluster.json file.")
    parser.add_argument("output_filtered_selected_csv", type=Path, help="Path to save the CSV with all filtered chains and selection marks.")
    parser.add_argument("output_selected_clusters_json", type=Path, help="Path to save the new cluster.json containing only selected representative clusters for resplitting.")
    parser.add_argument("--min_valid_proportion", type=float, default=0.7, help="Minimum seq_valid_proportion to keep a chain (default: 0.7).")
    parser.add_argument("--random_seed", type=int, default=66, help="Random seed for tie-breaking in representative selection (default: 66).")
    parser.add_argument("--log_file", type=Path, default=None, help="Optional path to a log file.")
    parser.add_argument("--log_level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="Logging level.")

    args = parser.parse_args()

    setup_logging(log_level_str=args.log_level, log_filepath=args.log_file)

    # --- 1. Load analysis results ---
    logger.info(f"Loading RNA chain analysis results from: {args.analysis_csv}")
    try:
        df_analysis = pd.read_csv(args.analysis_csv, dtype={'Resolution': str}) # Keep resolution as string initially
    except FileNotFoundError:
        logger.error(f"Analysis CSV file not found: {args.analysis_csv}")
        sys.exit(1)
    logger.info(f"Loaded {len(df_analysis)} chains from analysis CSV.")

    # --- 2. Filter by Pairing_Type == 'intramolecular' ---
    logger.info("Filtering for 'intramolecular' pairing type...")
    df_intra = df_analysis[df_analysis['Pairing_Type'] == 'intramolecular'].copy() # Use .copy() to avoid SettingWithCopyWarning
    logger.info(f"Found {len(df_intra)} chains with intramolecular pairing.")
    if df_intra.empty:
        logger.warning("No intramolecular chains found. Output files will be empty or reflect this.")
        # Create empty files and exit gracefully or proceed if downstream steps can handle empty inputs
        pd.DataFrame(columns=df_analysis.columns.tolist() + ['seq_valid_proportion', 'full_id', 'is_selected']).to_csv(args.output_filtered_selected_csv, index=False)
        with open(args.output_selected_clusters_json, 'w') as f: json.dump({}, f)
        logger.info("Empty output files created as no intramolecular chains were found.")
        sys.exit(0)


    # --- 3. Calculate seq_valid_proportion ---
    # This was step 2 in the notebook, applied after intramolecular filter
    df_intra = calculate_seq_valid_proportion(df_intra)

    # --- 4. Filter by seq_valid_proportion ---
    logger.info(f"Filtering by seq_valid_proportion > {args.min_valid_proportion}...")
    df_valid_len_intra = df_intra[df_intra['seq_valid_proportion'] > args.min_valid_proportion].copy()
    logger.info(f"Found {len(df_valid_len_intra)} chains after length proportion filter.")
    if df_valid_len_intra.empty:
        logger.warning("No chains remaining after seq_valid_proportion filter. Output files will be empty or reflect this.")
        pd.DataFrame(columns=df_intra.columns.tolist() + ['full_id', 'is_selected']).to_csv(args.output_filtered_selected_csv, index=False)
        with open(args.output_selected_clusters_json, 'w') as f: json.dump({}, f)
        logger.info("Empty output files created as no chains remained after proportion filter.")
        sys.exit(0)

    # --- 5. Create 'full_id' (PDB_ID + '_' + Chain_ID) ---
    df_valid_len_intra['full_id'] = df_valid_len_intra['PDB_ID'].astype(str) + '_' + df_valid_len_intra['Chain_ID'].astype(str)

    # --- 6. Load RNA3DB cluster.json and merge metadata ---
    df_cluster_json_meta = rna3db_cluster_json_to_df(args.cluster_json)
    # Rename columns from cluster.json to avoid clashes, e.g., if 'resolution' also exists there
    df_cluster_json_meta = df_cluster_json_meta.rename(columns={
        'sequence': 'rna3db_cluster_sequence', 
        'length': 'rna3db_cluster_length',
        'resolution': 'rna3db_cluster_resolution' # Keep original resolution from analysis for selection
    })
    
    # Merge. We need 'repr_pdb_chain' and metadata like 'rna3db_cluster_resolution' from df_cluster_json_meta
    # The 'full_id' from df_valid_len_intra should match 'pdb_chain_id_from_cluster'
    df_merged = pd.merge(
        df_valid_len_intra, 
        df_cluster_json_meta, 
        how='inner', 
        left_on='full_id', 
        right_on='pdb_chain_id_from_cluster'
    )
    logger.info(f"Merged analysis data with cluster.json metadata. Resulting chains: {len(df_merged)}.")
    if df_merged.empty:
        logger.warning("No chains matched between analysis results and cluster.json after filtering. Output files will be empty.")
        # Create empty files
        df_valid_len_intra['is_selected'] = 0 # Add the column for schema consistency
        df_valid_len_intra.to_csv(args.output_filtered_selected_csv, index=False)
        with open(args.output_selected_clusters_json, 'w') as f: json.dump({}, f)
        logger.info("Empty output files created as no chains matched during merge.")
        sys.exit(0)

    # The 'resolution' column from the original analysis (df_valid_len_intra) should be used for selection.
    # If df_cluster_json_meta also had a 'resolution' (now 'rna3db_cluster_resolution'), we keep the one from analysis.
    # The merge might create _x, _y suffixes if column names were identical and not renamed.
    # For selection, we will use the 'Resolution' column from the analysis CSV.

    # --- 7. Select representative chain for each 'repr_pdb_chain' group ---
    logger.info("Selecting representative chain for each 'repr_pdb_chain' group...")
    # Apply the selection function to each group
    # The result of apply here will be a DataFrame of the selected rows
    selected_rows_list = []
    for name, group in df_merged.groupby('repr_pdb_chain'):
        selected_chain_series = select_representative_chain(group.copy(), random_seed=args.random_seed) # Pass copy to avoid modifying group in place
        if not selected_chain_series.empty:
            selected_rows_list.append(selected_chain_series)
    
    if not selected_rows_list:
        logger.warning("No representative chains were selected. This is unexpected if previous steps had data.")
        df_merged['is_selected'] = 0
        df_merged.to_csv(args.output_filtered_selected_csv, index=False)
        with open(args.output_selected_clusters_json, 'w') as f: json.dump({}, f)
        sys.exit(0)

    df_selected_representatives = pd.DataFrame(selected_rows_list)
    
    # Add 'is_selected' column to the merged DataFrame
    df_merged['is_selected'] = 0
    # Mark the selected rows. Need to use index from df_selected_representatives to locate in df_merged.
    # The 'select_representative_chain' returns a Series which still has its original index from df_merged.
    df_merged.loc[df_selected_representatives.index, 'is_selected'] = 1
    
    logger.info(f"Marked {df_merged['is_selected'].sum()} chains as selected representatives.")

    # --- 8. Save the CSV with all filtered chains and selection marks ---
    # This CSV contains all chains that passed initial filters, plus the 'is_selected' mark.
    # It also includes merged metadata from cluster.json.
    # Drop the temporary 'resolution_numeric' column if it exists
    if 'resolution_numeric' in df_merged.columns:
        df_merged = df_merged.drop(columns=['resolution_numeric'])
    if 'pdb_chain_id_from_cluster' in df_merged.columns and 'full_id' in df_merged.columns: # pdb_chain_id_from_cluster is redundant with full_id after merge
        if 'pdb_chain_id_from_cluster' != 'full_id': # if they are not the same column already
             df_merged = df_merged.drop(columns=['pdb_chain_id_from_cluster'])


    df_merged.to_csv(args.output_filtered_selected_csv, index=False)
    logger.info(f"Filtered and selection-marked data saved to: {args.output_filtered_selected_csv}")

    # --- 9. Prepare and save the new cluster.json for resplitting ---
    # This new JSON will only contain the representative PDB clusters ('repr_pdb_chain')
    # for which a chain was actually selected.
    # The structure should be identical to the original cluster.json, but subsetted.
    
    # Get the set of 'repr_pdb_chain' for which we have a selected representative
    unique_selected_repr_ids = set(df_merged[df_merged['is_selected'] == 1]['repr_pdb_chain'])
    
    filter_rna3db_json_by_repr_ids(
        input_json_path=args.cluster_json,
        selected_repr_pdb_ids=unique_selected_repr_ids,
        output_json_path=args.output_selected_clusters_json
    )
    
    logger.info("Script finished successfully.")

if __name__ == "__main__":
    main()