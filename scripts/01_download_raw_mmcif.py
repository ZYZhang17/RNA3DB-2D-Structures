#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to download mmCIF files for PDB IDs extracted from an RNA3DB cluster.json file.
"""

import json
import os
import subprocess
import sys
from pathlib import Path
import argparse

def extract_pdb_ids_from_cluster_json(data):
    """
    Extracts all unique PDB IDs from the RNA3DB cluster.json structure.
    从 RNA3DB cluster.json 结构中提取所有唯一的 PDB ID。

    Args:
        data (dict): Parsed cluster.json data.
                     cluster.json 的解析数据。

    Returns:
        list: Sorted list of unique PDB IDs (lowercase).
              排序后的唯一 PDB ID 列表（小写）。
    """
    pdb_ids = set()
    # cluster.json structure: {component_id: {repr_pdb_id: {actual_pdb_chain_id: metadata}}}
    for component in data.values(): # Iterate over components
        for repr_pdb_cluster in component.values(): # Iterate over representative PDB clusters within a component
            for chain_id_key in repr_pdb_cluster.keys(): # Iterate over actual PDB chain IDs in the cluster
                # chain_id_key is like "1a2b_A"
                pdb_id = chain_id_key.split('_')[0].lower()
                pdb_ids.add(pdb_id)
    return sorted(list(pdb_ids))

def count_chain_ids_from_cluster_json(data):
    """
    Counts the total number of PDB chain IDs in the cluster.json structure.
    统计 cluster.json 结构中 PDB chain ID 的总数。

    Args:
        data (dict): Parsed cluster.json data.
                     cluster.json 的解析数据。
    Returns:
        int: Total count of chain IDs.
             Chain ID 的总数。
    """
    chain_count = 0
    for component in data.values():
        for repr_pdb_cluster in component.values():
            chain_count += len(repr_pdb_cluster.keys())
    return chain_count

def prepare_download_list(pdb_ids, output_file):
    """
    Prepares the list of PDB IDs in the format required for batch download.
    将 PDB ID 列表准备为批量下载所需的格式。

    Args:
        pdb_ids (list): List of PDB IDs.
                        PDB ID 列表。
        output_file (str or Path): Path to save the download list.
                                   下载列表的保存路径。
    """
    with open(output_file, 'w') as f:
        f.write(','.join(pdb_ids))
    print(f"Download list prepared: {output_file} with {len(pdb_ids)} IDs.")

def download_pdbs(download_list_file, output_dir, batch_download_script):
    """
    Downloads PDB files using the batch_download.sh script.
    使用 batch_download.sh 脚本下载 PDB 文件。

    Args:
        download_list_file (str or Path): Path to the file containing comma-separated PDB IDs.
                                          包含逗号分隔的 PDB ID 的文件路径。
        output_dir (str or Path): Directory to save downloaded mmCIF files.
                                  保存下载的 mmCIF 文件的目录。
        batch_download_script (str or Path): Path to the batch_download.sh script.
                                             batch_download.sh 脚本的路径。
    """
    if not os.path.exists(download_list_file) or os.path.getsize(download_list_file) == 0:
        print(f"Skipping download: PDB ID list file is empty or not found - {download_list_file}")
        return

    print(f"Downloading mmCIF files to {output_dir} using {batch_download_script}...")
    # Ensure the batch_download_script is executable
    # subprocess.run(["chmod", "+x", str(batch_download_script)], check=True) # Might not be needed if already executable
    
    # Command: bash BATCH_DOWNLOAD_SCRIPT -f download_list -o output_dir -c
    # The -c flag downloads .cif.gz files
    cmd = [
        "bash",
        str(batch_download_script),
        "-f",
        str(download_list_file),
        "-o",
        str(output_dir),
        "-c"  # Download mmCIF compressed files
    ]
    try:
        process = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("Download script stdout:")
        print(process.stdout)
        if process.stderr:
            print("Download script stderr:")
            print(process.stderr)
        print("Download command executed.")
    except subprocess.CalledProcessError as e:
        print(f"Error during PDB download:")
        print(f"Command: {' '.join(e.cmd)}")
        print(f"Return code: {e.returncode}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        sys.exit(1)


def uncompress_files(directory):
    """
    Uncompresses all .gz files in the specified directory.
    解压指定目录中的所有 .gz 文件。

    Args:
        directory (str or Path): Directory containing .gz files.
                                 包含 .gz 文件的目录。
    """
    print(f"Uncompressing files in {directory}...")
    found_gz = False
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.cif.gz'): # Specifically target .cif.gz
                found_gz = True
                file_path = Path(root) / file
                uncompressed_file_path = file_path.with_suffix('') # Removes .gz
                
                # Check if uncompressed file already exists and is newer or same size
                if uncompressed_file_path.exists():
                    # A simple check, could be more robust (e.g. check if gz is newer than uncompressed)
                    # For now, if uncompressed exists, assume it's fine or will be overwritten by gunzip -f
                    print(f"Uncompressed file {uncompressed_file_path} already exists. Gunzip will attempt to overwrite if older.")

                print(f"Uncompressing {file_path} to {uncompressed_file_path}")
                try:
                    # Using gunzip -f to force overwrite if the output file exists and is older
                    # and to suppress some warnings if the .gz file is corrupted (though it will still fail)
                    subprocess.run(["gunzip", "-f", str(file_path)], check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error uncompressing {file_path}: {e}")
                    # Decide if you want to exit or continue
                    # sys.exit(1)
                except FileNotFoundError:
                    print(f"Error: gunzip command not found. Please ensure gunzip is installed and in your PATH.")
                    sys.exit(1)
    if not found_gz:
        print("No .cif.gz files found to uncompress.")
    else:
        print("Uncompression attempt complete.")


def main():
    parser = argparse.ArgumentParser(description="Download mmCIF files based on an RNA3DB cluster.json.")
    parser.add_argument("cluster_json_file", type=Path, help="Path to the RNA3DB cluster.json file.")
    parser.add_argument("output_dir", type=Path, help="Directory to save downloaded and uncompressed mmCIF files.")
    parser.add_argument("--batch_script", type=Path, required=True, help="Path to the batch_download.sh script.")
    
    args = parser.parse_args()

    # Validate paths
    if not args.cluster_json_file.is_file():
        print(f"Error: Cluster JSON file not found: {args.cluster_json_file}")
        sys.exit(1)
    
    if not args.batch_script.is_file():
        print(f"Error: Batch download script not found: {args.batch_script}")
        sys.exit(1)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load cluster.json data
    print(f"Loading data from {args.cluster_json_file}...")
    with open(args.cluster_json_file, 'r') as f:
        cluster_data = json.load(f)
    
    # Extract PDB IDs
    pdb_ids_to_download = extract_pdb_ids_from_cluster_json(cluster_data)
    total_chain_ids = count_chain_ids_from_cluster_json(cluster_data)
    
    print(f"Total PDB chain IDs in cluster.json: {total_chain_ids}")
    print(f"Total unique PDB IDs to download: {len(pdb_ids_to_download)}")
    
    if not pdb_ids_to_download:
        print("No PDB IDs found to download. Exiting.")
        sys.exit(0)
        
    # Prepare download list
    # Place download list in a temporary location or within the output_dir
    download_list_path = args.output_dir / "pdb_download_list.txt"
    prepare_download_list(pdb_ids_to_download, download_list_path)
    
    # Download PDBs
    # The batch_download.sh script will download them as .cif.gz
    download_pdbs(download_list_path, args.output_dir, args.batch_script)
    
    # Uncompress files
    uncompress_files(args.output_dir)
    
    print("Download and uncompression process complete!")

if __name__ == "__main__":
    main()