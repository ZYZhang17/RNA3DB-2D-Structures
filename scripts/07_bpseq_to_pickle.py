#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to convert BPSEQ files into pickle format for RNA secondary structure prediction models.
This script processes BPSEQ files from a specified directory and creates a pickle file containing
RNA secondary structure data in the format required by UFold/ERNIE-RNA models.

此脚本将BPSEQ文件转换为RNA二级结构预测模型使用的pickle格式。
它处理指定目录中的BPSEQ文件，并创建一个包含RNA二级结构数据的pickle文件，
格式符合UFold/ERNIE-RNA模型的要求。
"""

import numpy as np
import os
import subprocess
import collections
import pickle
import random
import sys
import json
import logging
import argparse
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Any
from tqdm import tqdm

# --- Global Logger ---
logger = logging.getLogger(__name__)

# --- Define RNA_SS_data named tuple ---
RNA_SS_data = collections.namedtuple('RNA_SS_data', 'seq ss_label length name pairs')

# --- Logging Setup ---
def setup_logging(log_level_str: str = "INFO", log_filepath: Optional[Path] = None):
    """
    设置日志记录。
    Set up logging.
    
    Args:
        log_level_str (str): 日志级别字符串 (DEBUG, INFO, WARNING, ERROR, CRITICAL)
                             Log level string (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_filepath (Optional[Path]): 日志文件路径，如果为None则只输出到控制台
                                      Log file path, if None, output to console only
    """
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


def one_hot(seq: str) -> np.ndarray:
    """
    将RNA序列转换为one-hot编码。
    Convert RNA sequence to one-hot encoding.
    
    Args:
        seq (str): RNA序列，由'A', 'U', 'C', 'G'组成
                  RNA sequence composed of 'A', 'U', 'C', 'G'
                  
    Returns:
        np.ndarray: one-hot编码的序列，形状为(len(seq), 4)
                    One-hot encoded sequence with shape (len(seq), 4)
    """
    RNN_seq = seq
    BASES = 'AUCG'
    bases = np.array([base for base in BASES])
    feat = np.concatenate(
            [[(bases == base.upper()).astype(int)] if str(base).upper() in BASES else np.array([[0] * len(BASES)]) for base
            in RNN_seq])

    return feat


def clean_pair(pair_list: List[List[int]], seq: str) -> List[List[int]]:
    """
    清理配对列表，只保留有效的碱基对。
    Clean the pairing list, keeping only valid base pairs.
    
    Args:
        pair_list (List[List[int]]): 碱基对列表
                                     List of base pairs
        seq (str): RNA序列
                   RNA sequence
                   
    Returns:
        List[List[int]]: 清理后的碱基对列表
                         Cleaned list of base pairs
    """
    for item in pair_list:
        if seq[item[0]] == 'A' and seq[item[1]] == 'U':
            continue
        elif seq[item[0]] == 'C' and seq[item[1]] == 'G':
            continue
        elif seq[item[0]] == 'U' and seq[item[1]] == 'A':
            continue
        elif seq[item[0]] == 'G' and seq[item[1]] == 'C':
            continue
        else:
            logger.debug(f'Invalid base pair: {seq[item[0]]}+{seq[item[1]]}')
            pair_list.remove(item)
    return pair_list


def process_bpseq_files(
    bpseq_dir: Path, 
    output_pickle: Path, 
    split_json_path: Optional[Path] = None,
    dataset_name: str = "test_set",
    random_seed: int = 0,
    max_length: int = 600
) -> None:
    """
    处理BPSEQ文件并生成pickle格式的数据集。
    Process BPSEQ files and generate a pickle format dataset.
    
    Args:
        bpseq_dir (Path): BPSEQ文件所在目录
                          Directory containing BPSEQ files
        output_pickle (Path): 输出pickle文件的路径
                             Path to output pickle file
        split_json_path (Optional[Path]): 包含数据集划分信息的JSON文件路径
                                         Path to JSON file containing dataset split information
        dataset_name (str): 数据集名称，用于在split JSON中查找
                           Dataset name, used to look up in split JSON
        random_seed (int): 随机种子，用于打乱文件顺序
                          Random seed for shuffling files
        max_length (int): 序列的最大长度，超过此长度的序列将被跳过
                         Maximum sequence length, sequences longer than this will be skipped
    """
    logger.info(f"Processing BPSEQ files from directory: {bpseq_dir}")
    
    # 加载split JSON（如果提供）
    split_data = None
    if split_json_path:
        try:
            logger.info(f"Loading split information from: {split_json_path}")
            with open(split_json_path, 'r') as f:
                split_data = json.load(f)
            logger.info(f"Successfully loaded split data with keys: {list(split_data.keys())}")
        except Exception as e:
            logger.warning(f"Failed to load split JSON file: {e}")
            split_data = None
    
    # 获取所有BPSEQ文件
    all_files_1 = os.listdir(bpseq_dir)
    all_files = []
    for file in all_files_1:
        if not file.startswith('.'):
            all_files.append(file)
    
    # 随机打乱文件顺序
    logger.info(f"Found {len(all_files)} BPSEQ files. Shuffling with random seed {random_seed}.")
    random.seed(random_seed)
    random.shuffle(all_files)
    
    # 处理所有文件
    all_files_list = []
    skipped_files = 0
    error_files = 0
    
    logger.info(f"Starting to process BPSEQ files...")
    for index, item_file in enumerate(tqdm(all_files, desc="Processing BPSEQ files")):
        file_path = os.path.join(bpseq_dir, item_file)
        
        # 读取序列（第2列）
        t0 = subprocess.getstatusoutput(f'awk \'{{print $2}}\' {file_path}')
        if t0[0] != 0:
            logger.warning(f"Failed to extract sequence from {item_file}")
            error_files += 1
            continue
            
        seq = ''.join(t0[1].split('\n'))
        
        # 验证序列 - 与原始代码保持一致的验证逻辑
        error = False
        if len(seq) < max_length + 1:
            # 尝试从split JSON中获取和验证序列
            if split_data and dataset_name in split_data:
                for components in split_data[dataset_name].keys():
                    for cluster_id in split_data[dataset_name][components].keys():
                        for pdb_id in split_data[dataset_name][components][cluster_id].keys():
                            if pdb_id == item_file.split('.')[0]:
                                seq_json = split_data[dataset_name][components][cluster_id][pdb_id]['sequence']
                                if len(seq_json) == len(seq):
                                    allowed = {'A', 'G', 'C', 'U'}
                                    if any(char.upper() not in allowed for char in seq_json):
                                        error = True
                                    else: 
                                        seq = seq_json
                                else:
                                    error = True
        
        if error:
            logger.warning(f"Sequence validation failed for {item_file}")
            skipped_files += 1
            continue
            
        # 验证序列只包含有效碱基
        allowed = {'A', 'G', 'C', 'U'}
        if any(char.upper() not in allowed for char in seq):
            logger.warning(f"Invalid bases in sequence: {item_file}")
            skipped_files += 1
            continue
        
        # 生成one-hot编码
        try:
            one_hot_matrix = one_hot(seq.upper())
        except Exception as e:
            logger.error(f"Failed to create one-hot encoding for {item_file}: {e}")
            error_files += 1
            continue
        
        # 读取位置（第1列）和配对信息（第3列）
        t1 = subprocess.getstatusoutput(f'awk \'{{print $1}}\' {file_path}')
        t2 = subprocess.getstatusoutput(f'awk \'{{print $3}}\' {file_path}')
        
        if t1[0] == 0 and t2[0] == 0:
            pair_dict_all_list = [[int(item_tmp)-1, int(t2[1].split('\n')[index_tmp])-1] 
                                 for index_tmp, item_tmp in enumerate(t1[1].split('\n')) 
                                 if int(t2[1].split('\n')[index_tmp]) != 0]
        else:
            logger.warning(f"Failed to extract pairing information from {item_file}")
            pair_dict_all_list = []
        
        # 获取序列名称和长度
        seq_name = item_file
        seq_len = len(seq)
        
        # 创建配对字典（只保留left < right的配对）
        pair_dict_all = dict([item for item in pair_dict_all_list if item[0] < item[1]])
        
        # 处理符合长度要求的序列
        if seq_len > 1 and seq_len <= max_length:
            # 创建二级结构标签
            ss_label = np.zeros((seq_len, 3), dtype=int)
            ss_label[[*pair_dict_all.keys()], ] = [0, 1, 0]  # 左配对
            ss_label[[*pair_dict_all.values()], ] = [0, 0, 1]  # 右配对
            ss_label[np.where(np.sum(ss_label, axis=1) <= 0)[0], ] = [1, 0, 0]  # 未配对
            
            # 标准化为max_length长度
            one_hot_matrix_padded = np.zeros((max_length, 4))
            one_hot_matrix_padded[:seq_len, ] = one_hot_matrix
            
            ss_label_padded = np.zeros((max_length, 3), dtype=int)
            ss_label_padded[:seq_len, ] = ss_label
            ss_label_padded[np.where(np.sum(ss_label_padded, axis=1) <= 0)[0], ] = [1, 0, 0]
            
            # 创建样本
            sample_tmp = RNA_SS_data(
                seq=one_hot_matrix_padded,
                ss_label=ss_label_padded,
                length=seq_len,
                name=seq_name,
                pairs=pair_dict_all_list
            )
            all_files_list.append(sample_tmp)
        else:
            logger.debug(f"Skipping {item_file}: sequence length {seq_len} not in range (1, {max_length}]")
            skipped_files += 1
    
    # 保存为pickle文件
    logger.info(f"Total processed files: {len(all_files)}")
    logger.info(f"Valid samples: {len(all_files_list)}")
    logger.info(f"Skipped files: {skipped_files}")
    logger.info(f"Error files: {error_files}")
    
    # 确保输出目录存在
    output_pickle.parent.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Saving pickle file to: {output_pickle}")
    with open(output_pickle, "wb") as f:
        pickle.dump(all_files_list, f)
    
    logger.info(f"Processing completed. Pickle file saved to {output_pickle}")


def main():
    parser = argparse.ArgumentParser(description="Convert BPSEQ files to pickle format for RNA secondary structure prediction models.")
    parser.add_argument("--bpseq_dir", type=Path, 
                        default=Path("./data/analysis_output/06_bpseq/test_set"),
                        help="Directory containing BPSEQ files to process.")
    parser.add_argument("--output_pickle", type=Path, 
                        default=Path("./data/analysis_output/pickle/test.pickle"),
                        help="Path to save the output pickle file.")
    parser.add_argument("--split_json", type=Path, 
                        default=None,
                        help="Optional path to a JSON file containing dataset split information.")
    parser.add_argument("--dataset_name", type=str, 
                        default="test_set",
                        help="Dataset name to use when looking up in split JSON file.")
    parser.add_argument("--random_seed", type=int, 
                        default=66,
                        help="Random seed for shuffling files.")
    parser.add_argument("--max_length", type=int, 
                        default=600,
                        help="Maximum sequence length to process. Sequences longer than this will be skipped.")
    parser.add_argument("--log_file", type=Path, 
                        default=Path("./data/analysis_output/07_bpseq_to_pickle.log"),
                        help="Path to save the log file.")
    parser.add_argument("--log_level", type=str, 
                        default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level.")
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(log_level_str=args.log_level, log_filepath=args.log_file)
    
    # Ensure output directory exists
    if args.log_file:
        args.log_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Process BPSEQ files
    process_bpseq_files(
        bpseq_dir=args.bpseq_dir,
        output_pickle=args.output_pickle,
        split_json_path=args.split_json,
        dataset_name=args.dataset_name,
        random_seed=args.random_seed,
        max_length=args.max_length
    )


if __name__ == "__main__":
    main()