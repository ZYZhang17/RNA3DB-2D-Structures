#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to convert selected RNA sequences and their corrected dot-bracket structures
into BPSEQ format files, organized by dataset (train/valid/test).

This script performs operations similar to steps 11-12 in the original
RDL1.6.7-1.5_valid_chain_selection.ipynb notebook.
"""

import pandas as pd
import os
import argparse
import sys
from pathlib import Path
from tqdm import tqdm
import logging
from typing import Optional # For type hinting

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


def dotbracket_to_bpseq_string(sequence: str, structure: str, full_id_for_log: str) -> str:
    """
    Converts a sequence and its dot-bracket secondary structure into BPSEQ format string.
    Handles only '(' and ')' for pairing. Other characters in structure are treated as unpaired.
    Returns an empty string if sequence and structure lengths do not match.

    将序列及其点括号二级结构转换为 BPSEQ 格式字符串。
    仅处理 '(' 和 ')' 进行配对。结构中的其他字符被视为空配对。
    如果序列和结构长度不匹配，则返回空字符串。

    Args:
        sequence (str): The RNA sequence (e.g., 'Sequence_Canonical').
                        RNA 序列 (例如, 'Sequence_Canonical')。
        structure (str): The dot-bracket structure (e.g., 'Structure_Corrected').
                         点括号结构 (例如, 'Structure_Corrected')。
        full_id_for_log (str): The full ID of the chain for logging (e.g., PDBID_ChainID).
                               用于日志记录的链的完整 ID (例如, PDBID_ChainID)。
    Returns:
        str: BPSEQ formatted string, or empty string on error.
             BPSEQ 格式的字符串，如果出错则为空字符串。
    """
    if not isinstance(sequence, str) or not isinstance(structure, str):
        logger.warning(f"Chain {full_id_for_log}: Sequence or structure is not a string. Seq type: {type(sequence)}, Struct type: {type(structure)}. Skipping BPSEQ conversion for this chain.")
        return ""
        
    if len(sequence) != len(structure):
        logger.warning(f"Chain {full_id_for_log}: Sequence length ({len(sequence)}) and structure length ({len(structure)}) mismatch. Skipping BPSEQ conversion.")
        return ""

    bpseq_lines = []
    pairing_partner = [0] * len(sequence) # 0 means unpaired
    
    # Use a stack to find base pairs for parentheses
    # This simple stack-based approach works for non-pseudoknotted structures.
    # For pseudoknots (e.g. involving '[', '{'), a more complex algorithm would be needed.
    # This script currently only processes '()' pairs.
    stack = []
    for i, char_struct in enumerate(structure):
        if char_struct == '(':
            stack.append(i)
        elif char_struct == ')':
            if stack:
                j = stack.pop()
                # Store 1-based indices for BPSEQ
                pairing_partner[i] = j + 1
                pairing_partner[j] = i + 1
            else:
                # Unmatched closing parenthesis - treat as unpaired, log warning
                logger.warning(f"Chain {full_id_for_log}: Unmatched ')' at position {i+1} in structure. Treating as unpaired.")
                pairing_partner[i] = 0 
    
    if stack: # Unmatched opening parentheses
        for i in stack:
            logger.warning(f"Chain {full_id_for_log}: Unmatched '(' at position {i+1} in structure. Treating as unpaired.")
            pairing_partner[i] = 0

    # Generate BPSEQ lines
    for i in range(len(sequence)):
        # BPSEQ format: index base partner_index (1-based)
        bpseq_lines.append(f"{i+1} {sequence[i]} {pairing_partner[i]}")
        
    return "\n".join(bpseq_lines)


def process_dataframe_to_bpseq_files(
    df: pd.DataFrame, 
    base_output_dir: Path,
    seq_col: str = 'Sequence_Canonical', 
    struct_col: str = 'Structure_Corrected',
    id_col: str = 'full_id',
    dataset_col: str = 'dataset'
    ):
    """
    Iterates through a DataFrame, converts specified sequence and structure columns
    to BPSEQ format, and saves them into subdirectories named after the 'dataset' column.

    遍历 DataFrame，将指定的序列和结构列转换为 BPSEQ 格式，
    并将其保存到以 'dataset' 列命名的子目录中。

    Args:
        df (pd.DataFrame): DataFrame containing the RNA data.
                           包含 RNA 数据的 DataFrame。
        base_output_dir (Path): The base directory to save BPSEQ files.
                                保存 BPSEQ 文件的基础目录。
        seq_col (str): Name of the column containing sequences.
                       包含序列的列名。
        struct_col (str): Name of the column containing dot-bracket structures.
                          包含点括号结构的列名。
        id_col (str): Name of the column containing the unique ID for filename (e.g., 'full_id').
                      包含用于文件名的唯一 ID 的列名 (例如, 'full_id')。
        dataset_col (str): Name of the column specifying the dataset (e.g., 'train_set').
                           指定数据集的列名 (例如, 'train_set')。
    """
    if not base_output_dir.exists():
        logger.info(f"Base output directory {base_output_dir} does not exist. Creating it.")
        base_output_dir.mkdir(parents=True, exist_ok=True)

    required_cols = [seq_col, struct_col, id_col, dataset_col]
    if not all(col in df.columns for col in required_cols):
        missing = [col for col in required_cols if col not in df.columns]
        logger.error(f"DataFrame is missing required columns for BPSEQ generation: {missing}. Aborting.")
        return

    logger.info(f"Processing {len(df)} chains for BPSEQ conversion...")
    
    chains_converted = 0
    chains_skipped_error = 0

    for _, row in tqdm(df.iterrows(), total=df.shape[0], desc="Generating BPSEQ files"):
        chain_id = row[id_col]
        dataset = row[dataset_col]
        sequence = row[seq_col]
        structure = row[struct_col]

        if pd.isna(chain_id) or pd.isna(dataset) or pd.isna(sequence) or pd.isna(structure):
            logger.warning(f"Chain with ID '{chain_id if not pd.isna(chain_id) else 'Unknown'}': Contains NaN in required fields. Skipping BPSEQ conversion.")
            chains_skipped_error += 1
            continue
        
        # Create dataset-specific subdirectory
        dataset_dir = base_output_dir / str(dataset)
        dataset_dir.mkdir(parents=True, exist_ok=True)
        
        bpseq_content = dotbracket_to_bpseq_string(sequence, structure, str(chain_id))
        
        if bpseq_content:
            output_filename = f"{chain_id}.bpseq"
            output_filepath = dataset_dir / output_filename
            try:
                with open(output_filepath, 'w', encoding='utf-8') as f:
                    f.write(bpseq_content)
                chains_converted += 1
            except IOError as e:
                logger.error(f"Error writing BPSEQ file {output_filepath} for chain {chain_id}: {e}")
                chains_skipped_error += 1
        else:
            # Error/warning already logged by dotbracket_to_bpseq_string
            chains_skipped_error += 1
            
    logger.info(f"BPSEQ file generation complete for this set.")
    logger.info(f"Successfully converted {chains_converted} chains to BPSEQ.")
    logger.info(f"Skipped or errored during conversion for {chains_skipped_error} chains.")


def main():
    parser = argparse.ArgumentParser(description="Generate BPSEQ files from RNA sequence and structure data.")
    parser.add_argument("--input_selected_seqs_csv", type=Path, 
                        default=Path("./data/analysis_output/05_final_selected_dataset.csv"),
                        help="Path to the CSV file containing selected sequences and their dataset assignments (output of 05_assign_dataset_info.py).")
    parser.add_argument("--output_bpseq_dir", type=Path, 
                        default=Path("./data/analysis_output/06_bpseq"),
                        help="Base directory to save the generated BPSEQ files (subdirectories per dataset will be created).")
    parser.add_argument("--generate_less_600nt", action="store_true", 
                        help="If set, also generate a separate set of BPSEQ files for sequences with canonical length < 600nt.")
    parser.add_argument("--less_600nt_subdir_name", type=str, 
                        default="bpseq_less600nt", 
                        help="Subdirectory name within output_bpseq_dir for sequences < 600nt, if --generate_less_600nt is set.")
    parser.add_argument("--log_file", type=Path, 
                        default=Path("./data/analysis_output/06_generate_bpseq.log"), 
                        help="Optional path to a log file.")
    parser.add_argument("--log_level", type=str, 
                        default="INFO", 
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], 
                        help="Logging level.")
    args = parser.parse_args()

    setup_logging(log_level_str=args.log_level, log_filepath=args.log_file)

    # --- 1. Load the final selected sequences CSV ---
    logger.info(f"Loading final selected sequences data from: {args.input_selected_seqs_csv}")
    try:
        df_selected_seqs = pd.read_csv(args.input_selected_seqs_csv)
    except FileNotFoundError:
        logger.error(f"Input CSV file not found: {args.input_selected_seqs_csv}")
        sys.exit(1)
    
    if df_selected_seqs.empty:
        logger.warning("Input CSV is empty. No BPSEQ files will be generated.")
        sys.exit(0)
    
    logger.info(f"Loaded {len(df_selected_seqs)} selected sequences for BPSEQ generation.")

    # --- 2. Process all selected sequences ---
    logger.info("--- Generating BPSEQ files for ALL selected sequences ---")
    process_dataframe_to_bpseq_files(
        df=df_selected_seqs,
        base_output_dir=args.output_bpseq_dir,
        seq_col='Sequence_Canonical', # As per notebook
        struct_col='Structure_Corrected', # As per notebook
        id_col='full_id', # As per notebook
        dataset_col='dataset' # As per notebook
    )

    # --- 3. Optionally, process sequences with length < 600nt ---
    if args.generate_less_600nt:
        logger.info("--- Generating BPSEQ files for selected sequences with canonical length < 600nt ---")
        
        # 使用与notebook一致的长度列 (rna3db_seq_length)
        if 'rna3db_seq_length' in df_selected_seqs.columns:
            # 优先使用rna3db_seq_length列（与RDL1.6.7-1.5 notebook一致）
            df_less_600nt = df_selected_seqs[df_selected_seqs['rna3db_seq_length'] < 600].copy()
            logger.info(f"Using 'rna3db_seq_length' column for length filtering. Found {len(df_less_600nt)} sequences < 600nt.")
        elif 'Sequence_Canonical' in df_selected_seqs.columns:
            # 备选方案：计算Sequence_Canonical长度（与原代码逻辑一致）
            logger.warning("Column 'rna3db_seq_length' not found. Using calculated length from 'Sequence_Canonical' as fallback.")
            df_selected_seqs['temp_canonical_length'] = df_selected_seqs['Sequence_Canonical'].apply(
                lambda x: len(x) if isinstance(x, str) else 0
            )
            df_less_600nt = df_selected_seqs[df_selected_seqs['temp_canonical_length'] < 600].copy()
            df_less_600nt.drop(columns=['temp_canonical_length'], inplace=True)
            logger.info(f"Found {len(df_less_600nt)} sequences with canonical length < 600nt.")
        else:
            logger.error("Neither 'rna3db_seq_length' nor 'Sequence_Canonical' column found for length filtering. Skipping <600nt processing.")
            df_less_600nt = pd.DataFrame()  # 空DataFrame，跳过后续处理

        if not df_less_600nt.empty:
            output_dir_less_600nt = args.output_bpseq_dir / args.less_600nt_subdir_name
            process_dataframe_to_bpseq_files(
                df=df_less_600nt,
                base_output_dir=output_dir_less_600nt,
                seq_col='Sequence_Canonical',
                struct_col='Structure_Corrected',
                id_col='full_id',
                dataset_col='dataset'
            )
        else:
            logger.info("No sequences found for <600nt processing or processing was skipped due to missing columns.")
    
    logger.info("Script finished successfully.")

if __name__ == "__main__":
    main()