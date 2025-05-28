#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to assign dataset information (train_set, valid_set, test_set) to filtered
and selected RNA chains from the previous step. This script reads the filtered and
selected chains CSV and the resplit JSON, then merges them to create final CSV files.

此脚本将数据集信息（train_set, valid_set, test_set）分配给上一步过滤和选择的RNA链。
脚本读取过滤和选择的链CSV以及resplit JSON，然后合并它们以创建最终的CSV文件。
"""

import pandas as pd
import json
import sys
import logging
from pathlib import Path
from typing import Optional

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


def parse_resplit_json_to_dataframe(resplit_json_path: Path) -> pd.DataFrame:
    """
    解析 resplit JSON 文件为 DataFrame，与 notebook 中的 read_dataset_partition 函数一致。
    Parse resplit JSON file into a DataFrame, consistent with the read_dataset_partition function in the notebook.
    
    Args:
        resplit_json_path (Path): resplit JSON 文件路径
                                  Path to the resplit JSON file.
        
    Returns:
        pd.DataFrame: 包含 dataset, component, repr_pdb_chain, pdb_chain_id 列的 DataFrame
                     DataFrame with dataset, component, repr_pdb_chain, pdb_chain_id columns.
    """
    logger.info(f"Parsing resplit JSON from: {resplit_json_path}")
    with open(resplit_json_path, 'r') as f:
        data = json.load(f)
    
    rows = []
    for dataset, components in data.items():
        for component, pdbs in components.items():
            for repr_pdb_chain, chains in pdbs.items():
                for pdb_chain_id, metas in chains.items():
                    row = {
                        "dataset": dataset,
                        "component": component,
                        "repr_pdb_chain": repr_pdb_chain,
                        "pdb_chain_id": pdb_chain_id
                    }
                    # 添加所有元数据 (Add all metadata)
                    row.update(metas)
                    rows.append(row)
    
    split_df = pd.DataFrame(rows)
    logger.info(f"Parsed resplit JSON into DataFrame with {len(split_df)} rows.")
    return split_df


def assign_dataset_info(filtered_selected_csv: Path, resplit_json: Path, 
                       output_final_all_info_csv: Path, output_final_selected_dataset_csv: Path):
    """
    分配数据集信息，使用与 notebook 相同的合并方式。
    Assign dataset information using the same merging method as in the notebook.
    
    Args:
        filtered_selected_csv (Path): 过滤和选择的链CSV文件路径
                                      Path to the filtered and selected chains CSV file.
        resplit_json (Path): resplit JSON 文件路径
                            Path to the resplit JSON file.
        output_final_all_info_csv (Path): 输出包含所有链的最终CSV文件路径
                                         Path to the output final CSV file with all chains.
        output_final_selected_dataset_csv (Path): 最终用于训练的RNA3DB dataset，输出仅包含选定序列的最终CSV文件路径
                                              Path to the output final CSV file with only selected sequences.
    """
    logger.info(f"Loading filtered and selected chains from: {filtered_selected_csv}")
    df = pd.read_csv(filtered_selected_csv)
    logger.info(f"Loaded {len(df)} chains from filtered CSV.")
    
    # 解析 resplit JSON 为 DataFrame（与 notebook 一致）
    # Parse resplit JSON to DataFrame (consistent with notebook)
    split_df = parse_resplit_json_to_dataframe(resplit_json)
    
    # 使用与 notebook 相同的合并方式：left_on='full_id', right_on='pdb_chain_id'
    # Use the same merging method as in the notebook: left_on='full_id', right_on='pdb_chain_id'
    logger.info("Merging dataset information using full_id and pdb_chain_id...")
    final_df = pd.merge(df, split_df, left_on='full_id', right_on='pdb_chain_id', how='left')
    
    # 检查未分配的链
    # Check for unassigned chains
    unassigned_count = final_df['dataset'].isna().sum()
    if unassigned_count > 0:
        logger.warning(f"Found {unassigned_count} chains without dataset assignment.")
        # 将未分配的标记为 'unassigned'
        # Mark unassigned chains as 'unassigned'
        final_df['dataset'] = final_df['dataset'].fillna('unassigned')
    
    logger.info(f"Dataset assignment complete. Final DataFrame has {len(final_df)} rows.")
    
    # 显示数据集分布
    # Display dataset distribution
    dataset_counts = final_df['dataset'].value_counts()
    logger.info(f"Dataset distribution: {dataset_counts.to_dict()}")
    
    # 保存包含所有信息的 CSV
    # Save CSV with all information
    logger.info(f"Saving final all info CSV to: {output_final_all_info_csv}")
    final_df.to_csv(output_final_all_info_csv, index=False)
    
    # 保存只包含选定序列的 CSV（is_selected == 1）
    # Save CSV with only selected sequences (is_selected == 1)
    final_selected_df = final_df[final_df['is_selected'] == 1].copy()
    logger.info(f"Saving final selected sequences CSV to: {output_final_selected_dataset_csv}")
    logger.info(f"Selected sequences count: {len(final_selected_df)}")
    final_selected_df.to_csv(output_final_selected_dataset_csv, index=False)
    
    # 显示选定序列的数据集分布
    # Display dataset distribution for selected sequences
    selected_dataset_counts = final_selected_df['dataset'].value_counts()
    logger.info(f"Selected sequences dataset distribution: {selected_dataset_counts.to_dict()}")
    
    logger.info("Dataset assignment completed successfully.")


def main():
    """
    Main function to parse command line arguments and run the script.
    解析命令行参数并运行脚本的主函数。
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Assign dataset information to filtered and selected RNA chains.")
    parser.add_argument("--filtered_selected_csv", type=Path, 
                        default=Path("./data/analysis_output/03_filtered_selected_chains.csv"),
                        help="Path to the filtered and selected chains CSV file from 03_filter_select_chains.py.")
    parser.add_argument("--resplit_json", type=Path, 
                        default=Path("./data/analysis_output/04_resplit.json"),
                        help="Path to the resplit JSON file from 04_resplit_dataset.sh.")
    parser.add_argument("--output_final_all_info_csv", type=Path, 
                        default=Path("./data/analysis_output/05_final_all_info.csv"),
                        help="Path to save the CSV with all chains and dataset information.")
    parser.add_argument("--output_final_selected_dataset_csv", type=Path, 
                        default=Path("./data/analysis_output/05_final_selected_dataset.csv"),
                        help="Path to save the CSV with only selected sequences and dataset information.")
    parser.add_argument("--log_file", type=Path, default="./data/analysis_output/05_assgin_dataset_info.log", 
                        help="Optional path to a log file.")
    parser.add_argument("--log_level", type=str, default="INFO", 
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], 
                        help="Logging level.")
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(log_level_str=args.log_level, log_filepath=args.log_file)
    
    # Ensure output directories exist
    args.output_final_all_info_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_final_selected_dataset_csv.parent.mkdir(parents=True, exist_ok=True)
    
    # Run the main function
    assign_dataset_info(
        filtered_selected_csv=args.filtered_selected_csv,
        resplit_json=args.resplit_json,
        output_final_all_info_csv=args.output_final_all_info_csv,
        output_final_selected_dataset_csv=args.output_final_selected_dataset_csv
    )


if __name__ == "__main__":
    main()