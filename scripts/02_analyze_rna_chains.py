#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RNA 3D 结构分析脚本 (RDL1.6.7-1.4)
RNA 3D Structure Analysis Script (RDL1.6.7-1.4)

功能/Function:
    遍历指定目录下的所有 mmCIF (.cif) 文件，利用 rnapolis-py 库分析每个文件，
    重点关注 RNA 链，并生成包含权威序列、RNApolis 原始输出以及
    修正后序列/结构的 CSV 报告。
    **修正：基于最新的 RNApolis 映射逻辑进行序列/结构修正。**

    Traverse all mmCIF (.cif) files in a given directory, analyze each file using rnapolis-py,
    focus on RNA chains, and generate a CSV report containing authoritative sequences,
    raw RNApolis output, and corrected sequences/structures.

输入/Input:
    - 一个包含 mmCIF 文件的目录路径。
      Directory path containing mmCIF files.
    - rna3db 代码库路径 (通过 sys.path 添加)。
      Path to rna3db codebase (added to sys.path).
    - modifications_cache.json 文件路径。
      Path to modifications_cache.json.

输出/Output:
    - 一个 CSV 文件，包含以下列:
      A CSV file with the following columns:
        - PDB_ID
        - Chain_ID (权威链 ID, e.g., pdb_strand_id)
        - Sequence_Canonical (权威序列, 仅含 ACGUN)
        - Sequence_RNApolis_Raw (RNApolis 原始序列, find_gaps=False)
        - Structure_RNApolis_Raw (RNApolis 原始点括号结构)
        - Sequence_Corrected (修正后序列, 与权威等长, 含 -/ACGUN)
        - Structure_Corrected (修正后结构, 与权威等长, 含 .)
        - Resolution (分辨率)
        - Pairing_Type (分子内/分子间/无)
        - Pairing_Partners (配对链 ID 列表)
    - 一个日志文件。
      A log file.

依赖/Dependencies:
    - rnapolis-py
    - mmcif
    - tqdm
    - pandas
    - rna3db (需要能访问其代码)

用法/Usage:
    python 02_analyze_rna_chains.py \
        --input_dir <mmcif目录> \
        --output_csv <输出csv路径> \
        --mod_cache <modifications_cache.json路径> \
        --log_file <日志文件路径> \
        [--max_workers <并行进程数>]

    Example:
    python 02_analyze_rna_chains.py \
        --input_dir ./data/downloaded_mmcif \
        --output_csv ./data/analysis_output/02_rna_chain_analysis.csv \
        --mod_cache ./utils/modifications_cache.json \
        --log_file ./data/analysis_output/02_analyze_rna_chains.log \
        --max_workers 8
"""

import csv
import logging
import os
import sys
import time
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Set, Tuple, Optional, Any
from functools import partial

# ----------------------------------------------------
# 导入 rnapolis 组件 / Import rnapolis components
from rnapolis.parser import read_3d_structure
from rnapolis.annotator import extract_secondary_structure
from rnapolis.common import Structure2D, Molecule, ResidueAuth, ResidueLabel, BasePair, BpSeq
from rnapolis.tertiary import Residue3D, Mapping2D3D, Structure3D
from rnapolis.util import handle_input_file

# 导入 mmcif 解析库 / Import mmcif parser
from mmcif.io.IoAdapterPy import IoAdapterPy
from tqdm import tqdm
import pandas as pd

# --- 尝试导入 ModificationHandler ---
# Try to import ModificationHandler from utils
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

try:
    from utils.modification_utils import ModificationHandler
    print("成功导入 ModificationHandler / Successfully imported ModificationHandler")
    USING_REAL_MODIFICATION_HANDLER = True
except ImportError as e:
    print(f"错误：无法导入ModificationHandler: {e} / Error: Cannot import ModificationHandler: {e}")
    USING_REAL_MODIFICATION_HANDLER = False
# ------------------------------------

def setup_logging(log_level_str: str, log_filepath: str):
    """
    设置日志记录器（中英文日志）
    Setup logging to file and console (bilingual log).
    """
    log_level = getattr(logging, log_level_str.upper(), logging.INFO)
    log_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    if root_logger.hasHandlers():
        root_logger.handlers.clear()
    try:
        file_handler = logging.FileHandler(log_filepath, mode='w', encoding='utf-8')
        file_handler.setFormatter(log_formatter)
        file_handler.setLevel(log_level)
        root_logger.addHandler(file_handler)
    except Exception as e:
        print(f"错误：无法设置日志文件 handler 到 {log_filepath}: {e}\n"
              f"Error: Cannot set log file handler to {log_filepath}: {e}")
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)
    console_handler.setLevel(log_level)
    root_logger.addHandler(console_handler)

logger = logging.getLogger(__name__)

def _parse_pdbx_poly_seq_scheme_detailed(
    cif_file_path: str,
    modification_handler: 'ModificationHandler'
) -> Optional[Dict[str, List[Dict]]]:
    """
    解析 mmCIF 文件中的权威序列信息
    Parse authoritative sequence info from mmCIF file.

    Args:
        cif_file_path (str): mmCIF 文件路径 / Path to mmCIF file.
        modification_handler (ModificationHandler): 修饰符处理器 / Modification handler.

    Returns:
        dict or None: {chain_id: [residue_dict, ...]} 或 None
    """
    scheme_data_by_chain = defaultdict(list)
    try:
        io_adapter = IoAdapterPy()
        data_containers = io_adapter.readFile(cif_file_path)
        if not data_containers:
            logger.error(f"文件 {os.path.basename(cif_file_path)}: 无法读取 mmCIF 数据容器。/ "
                         f"File {os.path.basename(cif_file_path)}: Cannot read mmCIF data container.")
            return None
        data_block = data_containers[0]
        scheme_obj = data_block.getObj("pdbx_poly_seq_scheme")
        if not scheme_obj:
            logger.info(f"文件 {os.path.basename(cif_file_path)}: 未找到 'pdbx_poly_seq_scheme' 对象。/ "
                        f"File {os.path.basename(cif_file_path)}: 'pdbx_poly_seq_scheme' object not found.")
            return {}

        attributes = scheme_obj.getAttributeList()
        required_attrs_base = ["pdb_strand_id", "seq_id", "mon_id"]
        auth_seq_num_attr = None
        if "pdb_seq_num" in attributes: auth_seq_num_attr = "pdb_seq_num"
        elif "auth_seq_num" in attributes: auth_seq_num_attr = "auth_seq_num"
        else:
            logger.error(f"文件 {os.path.basename(cif_file_path)}: 解析 _pdbx_poly_seq_scheme 时缺少作者序列号属性。/ "
                         f"File {os.path.basename(cif_file_path)}: Missing author sequence number attribute in _pdbx_poly_seq_scheme.")
            return None

        required_attrs = required_attrs_base + [auth_seq_num_attr]
        if not all(attr in attributes for attr in required_attrs):
            missing = [attr for attr in required_attrs if attr not in attributes]
            logger.error(f"文件 {os.path.basename(cif_file_path)}: 解析 _pdbx_poly_seq_scheme 时缺少必要属性: {missing}. / "
                         f"File {os.path.basename(cif_file_path)}: Missing required attributes in _pdbx_poly_seq_scheme: {missing}.")
            return None

        chain_idx = attributes.index("pdb_strand_id")
        seq_id_idx = attributes.index("seq_id")
        mon_id_idx = attributes.index("mon_id")
        auth_seq_num_idx = attributes.index(auth_seq_num_attr)
        ins_code_idx = attributes.index("pdb_ins_code") if "pdb_ins_code" in attributes else -1

        for row in scheme_obj.getRowList():
            try:
                chain_id = row[chain_idx]
                seq_id_val = row[seq_id_idx]
                seq_id = int(seq_id_val) if seq_id_val not in ['.', '?'] else -1
                mon_id = row[mon_id_idx]
                auth_seq_num_val = row[auth_seq_num_idx]
                auth_seq_num = int(auth_seq_num_val) if auth_seq_num_val not in ['.', '?'] else -1
                ins_code_raw = row[ins_code_idx] if ins_code_idx != -1 else '?'
                ins_code = ins_code_raw if ins_code_raw not in ['?', '.'] else None
                std_one_letter_code = modification_handler.rna_letters_3to1(mon_id)
                scheme_data_by_chain[chain_id].append({
                    "seq_id": seq_id, "mon_id": mon_id, "std_one_letter_code": std_one_letter_code,
                    "auth_seq_num": auth_seq_num, "pdb_ins_code": ins_code
                })
            except (ValueError, IndexError) as e:
                logger.warning(f"文件 {os.path.basename(cif_file_path)}: 解析 _pdbx_poly_seq_scheme 行时出错: {row} - {e} / "
                               f"File {os.path.basename(cif_file_path)}: Error parsing _pdbx_poly_seq_scheme row: {row} - {e}")
                continue

        for chain_id_key in scheme_data_by_chain:
            scheme_data_by_chain[chain_id_key].sort(key=lambda x: x["seq_id"])
        return dict(scheme_data_by_chain)

    except FileNotFoundError:
        logger.error(f"文件 {os.path.basename(cif_file_path)}: 尝试解析权威序列时文件未找到。/ "
                     f"File {os.path.basename(cif_file_path)}: File not found when parsing authoritative sequence.")
        return None
    except Exception as e:
        logger.error(f"文件 {os.path.basename(cif_file_path)}: 解析 _pdbx_poly_seq_scheme 时发生意外错误: {e} / "
                     f"File {os.path.basename(cif_file_path)}: Unexpected error parsing _pdbx_poly_seq_scheme: {e}", exc_info=False)
        return None

def parse_multiline_dot_bracket(dot_bracket_multiline: str, pdb_id: str) -> Dict[str, Tuple[str, str]]:
    """
    解析多行点括号结构字符串
    Parse multiline dot-bracket structure string.

    Args:
        dot_bracket_multiline (str): 多行点括号字符串 / Multiline dot-bracket string.
        pdb_id (str): PDB ID.

    Returns:
        dict: {chain_id: (sequence, structure)}
    """
    parsed_data = {}
    if not dot_bracket_multiline or not dot_bracket_multiline.strip():
        logger.warning(f"文件 {pdb_id}: 传入的点括号字符串为空。/ File {pdb_id}: Input dot-bracket string is empty.")
        return parsed_data

    strands_data = dot_bracket_multiline.strip().split('\n>')
    if strands_data and strands_data[0].startswith('>'):
        strands_data[0] = strands_data[0][1:]

    for i, strand_block in enumerate(strands_data):
        if not strand_block.strip(): continue
        lines = strand_block.strip().split('\n')
        if len(lines) == 3:
            strand_id_raw = lines[0].strip()
            sequence = lines[1].strip()
            structure = lines[2].strip()
            strand_id_key = strand_id_raw[len("strand_"):] if strand_id_raw.lower().startswith("strand_") else strand_id_raw
            if strand_id_key in parsed_data:
                old_seq_len = len(parsed_data[strand_id_key][0])
                new_seq_len = len(sequence)
                if new_seq_len > old_seq_len:
                    logger.info(f"文件 {pdb_id}: 重复链ID '{strand_id_key}'，新数据更长，已覆盖。/ "
                                f"File {pdb_id}: Duplicate chain ID '{strand_id_key}', new data is longer, replaced.")
                    parsed_data[strand_id_key] = (sequence, structure)
                else:
                    logger.info(f"文件 {pdb_id}: 重复链ID '{strand_id_key}'，新数据不长于旧数据，保留旧数据。/ "
                                f"File {pdb_id}: Duplicate chain ID '{strand_id_key}', new data not longer, keep old.")
            else:
                parsed_data[strand_id_key] = (sequence, structure)
        elif len(lines) == 1 and not lines[0]: continue
        else:
            logger.warning(f"文件 {pdb_id}: 解析点括号块时遇到意外格式: {lines} / "
                           f"File {pdb_id}: Unexpected format when parsing dot-bracket block: {lines}")
    return parsed_data

def _parse_resolution(cif_file_path: str) -> str:
    """
    解析 mmCIF 文件中的分辨率
    Parse resolution from mmCIF file.

    Args:
        cif_file_path (str): mmCIF 文件路径 / Path to mmCIF file.

    Returns:
        str: 分辨率字符串 / Resolution string.
    """
    try:
        io_adapter = IoAdapterPy()
        data_containers = io_adapter.readFile(cif_file_path)
        if not data_containers: return "N/A"
        data_block = data_containers[0]
        refine_obj = data_block.getObj("refine")
        if not refine_obj: return "N/A"
        attributes = refine_obj.getAttributeList(); res_attr = "ls_d_res_high"
        if res_attr not in attributes: return "N/A"
        res_idx = attributes.index(res_attr); rows = refine_obj.getRowList()
        if not rows: return "N/A"
        resolution_val = rows[0][res_idx]
        return resolution_val if resolution_val not in ['?', '.'] else "N/A"
    except Exception:
        return "N/A"

def process_file(
    cif_path: str,
    modification_handler: 'ModificationHandler'
) -> List[List[str]]:
    """
    处理单个 mmCIF 文件，提取 RNA 链信息并生成结果行
    Process a single mmCIF file, extract RNA chain info, and generate result rows.

    Args:
        cif_path (str): mmCIF 文件路径 / Path to mmCIF file.
        modification_handler (ModificationHandler): 修饰符处理器 / Modification handler.

    Returns:
        list: 每条 RNA 链一行的结果 / List of result rows for each RNA chain.
    """
    file_rows = []
    pdb_id = os.path.splitext(os.path.basename(cif_path))[0]
    logger.info(f"开始处理文件: {pdb_id} / Start processing file: {pdb_id}")
    temp_cif_handle = None
    rna_chain_auth_ids = set()
    existing_3d_residues_by_auth_id: Dict[Tuple, Residue3D] = {}
    canonical_residues_by_pdb_strand_id: Optional[Dict[str, List[Dict]]] = None
    structure3d: Optional[Structure3D] = None
    structure2d_obj: Optional[Structure2D] = None
    mapping: Optional[Mapping2D3D] = None
    rnapolis_parsed_data = {}
    residue_to_rnapolis_idx = {}
    rnapolis_idx_to_chars = {}
    pairing_partners = defaultdict(set)
    rna_structure_type = {}
    resolution = "N/A"

    try:
        temp_cif_handle = handle_input_file(cif_path)
        temp_cif_path = temp_cif_handle.name

        canonical_residues_by_pdb_strand_id = _parse_pdbx_poly_seq_scheme_detailed(
            temp_cif_path, modification_handler
        )
        if canonical_residues_by_pdb_strand_id is None:
            logger.error(f"文件 {pdb_id}: 无法解析权威序列信息，跳过。/ "
                         f"File {pdb_id}: Cannot parse authoritative sequence info, skipped.")
            return file_rows

        resolution = _parse_resolution(temp_cif_path)

        temp_cif_handle.seek(0)
        structure3d = read_3d_structure(temp_cif_handle, model=None, nucleic_acid_only=False)
        if not structure3d or not structure3d.residues:
            logger.warning(f"文件 {pdb_id}: RNApolis 未找到3D残基。/ "
                           f"File {pdb_id}: RNApolis did not find 3D residues.")
            return file_rows

        chain_molecule_types: Dict[str, set] = defaultdict(set)
        for residue in structure3d.residues:
            auth_chain_id = residue.auth.chain if residue.auth else None
            if auth_chain_id:
                chain_molecule_types[auth_chain_id].add(residue.molecule_type)
                if residue.auth:
                    auth_id = (residue.auth.chain, residue.auth.number, residue.auth.icode)
                    existing_3d_residues_by_auth_id[auth_id] = residue
        for chain_id, types in chain_molecule_types.items():
            if Molecule.RNA in types: rna_chain_auth_ids.add(chain_id)

        if not rna_chain_auth_ids:
            logger.info(f"文件 {pdb_id}: 未找到 RNA 链。/ File {pdb_id}: No RNA chain found.")
            return file_rows

        try:
            result_tuple = extract_secondary_structure(structure3d, model=None, find_gaps=False, all_dot_brackets=False)
            if result_tuple and isinstance(result_tuple[0], Structure2D):
                structure2d_obj = result_tuple[0]
                base_pairs = structure2d_obj.baseInteractions.basePairs if structure2d_obj.baseInteractions else []
                stackings = structure2d_obj.baseInteractions.stackings if structure2d_obj.baseInteractions else []
                mapping = Mapping2D3D(structure3d, base_pairs, stackings, False)
                rnapolis_parsed_data = parse_multiline_dot_bracket(structure2d_obj.dotBracket, pdb_id)

                if mapping:
                    if hasattr(mapping, 'bpseq_index_to_residue_map') and mapping.bpseq_index_to_residue_map:
                        temp_map = mapping.bpseq_index_to_residue_map
                        residue_to_rnapolis_idx = {v: k for k, v in temp_map.items()}
                    else:
                        logger.warning(f"文件 {pdb_id}: mapping.bpseq_index_to_residue_map 不可用。/ "
                                       f"File {pdb_id}: mapping.bpseq_index_to_residue_map not available.")

                    if hasattr(mapping, 'bpseq') and isinstance(mapping.bpseq, BpSeq) and \
                       hasattr(mapping.bpseq, 'dot_bracket') and mapping.bpseq.dot_bracket and \
                       hasattr(mapping.bpseq.dot_bracket, 'sequence') and \
                       hasattr(mapping.bpseq.dot_bracket, 'structure'):
                        internal_seq = mapping.bpseq.dot_bracket.sequence
                        internal_struct = mapping.bpseq.dot_bracket.structure
                        if len(internal_seq) == len(internal_struct):
                            for i in range(len(internal_seq)): rnapolis_idx_to_chars[i + 1] = (internal_seq[i], internal_struct[i])
                        else:
                            logger.error(f"文件 {pdb_id}: mapping.bpseq.dot_bracket 序列/结构长度不匹配。/ "
                                         f"File {pdb_id}: mapping.bpseq.dot_bracket sequence/structure length mismatch.")
                    else:
                        logger.warning(f"文件 {pdb_id}: mapping.bpseq.dot_bracket 不可用，回退构建 rnapolis_idx_to_chars。/ "
                                       f"File {pdb_id}: mapping.bpseq.dot_bracket not available, fallback to build rnapolis_idx_to_chars.")
                        current_idx = 1
                        sorted_strand_keys = sorted(rnapolis_parsed_data.keys())
                        temp_s, temp_db = "", ""
                        for key in sorted_strand_keys: seq, struct = rnapolis_parsed_data[key]; temp_s += seq; temp_db += struct
                        if len(temp_s) == len(temp_db):
                            for i in range(len(temp_s)): rnapolis_idx_to_chars[current_idx] = (temp_s[i], temp_db[i]); current_idx += 1
                        else:
                            logger.error(f"文件 {pdb_id}: 回退构建 rnapolis_idx_to_chars 时长度不匹配。/ "
                                         f"File {pdb_id}: Fallback build rnapolis_idx_to_chars length mismatch.")
                else:
                    logger.warning(f"文件 {pdb_id}: Mapping2D3D 对象未创建。/ "
                                   f"File {pdb_id}: Mapping2D3D object not created.")
            else:
                logger.warning(f"文件 {pdb_id}: RNApolis 未能提取二级结构。/ "
                               f"File {pdb_id}: RNApolis failed to extract secondary structure.")
        except Exception as e_rnapolis:
            logger.error(f"文件 {pdb_id}: RNApolis 分析时出错: {e_rnapolis} / "
                         f"File {pdb_id}: Error during RNApolis analysis: {e_rnapolis}", exc_info=False)

        if structure2d_obj:
            rna_chain_has_pairs = defaultdict(bool)
            base_pairs_from_rnapolis: List[BasePair] = structure2d_obj.baseInteractions.basePairs
            for bp in base_pairs_from_rnapolis:
                res1_3d = structure3d.find_residue(bp.nt1.label, bp.nt1.auth)
                res2_3d = structure3d.find_residue(bp.nt2.label, bp.nt2.auth)
                if res1_3d and res1_3d.auth and res2_3d and res2_3d.auth:
                    nt1_auth_chain, nt2_auth_chain = res1_3d.auth.chain, res2_3d.auth.chain
                    if nt1_auth_chain and nt2_auth_chain:
                        pairing_partners[nt1_auth_chain].add(nt2_auth_chain)
                        pairing_partners[nt2_auth_chain].add(nt1_auth_chain)
                        if nt1_auth_chain in rna_chain_auth_ids:
                            rna_chain_has_pairs[nt1_auth_chain] = True
                            if nt1_auth_chain != nt2_auth_chain: rna_structure_type[nt1_auth_chain] = "intermolecular"
                        if nt2_auth_chain in rna_chain_auth_ids:
                            rna_chain_has_pairs[nt2_auth_chain] = True
                            if nt1_auth_chain != nt2_auth_chain and rna_structure_type.get(nt2_auth_chain) != "intermolecular":
                                rna_structure_type[nt2_auth_chain] = "intermolecular"
            for rna_auth_chain_id in rna_chain_auth_ids:
                 if rna_auth_chain_id not in rna_structure_type:
                     rna_structure_type[rna_auth_chain_id] = "intramolecular" if rna_chain_has_pairs[rna_auth_chain_id] else "no_secondary_structure_found"
        else:
             for rna_auth_chain_id in rna_chain_auth_ids: rna_structure_type[rna_auth_chain_id] = "N/A (RNApolis failed)"

        for pdb_strand_id, canonical_residue_list in canonical_residues_by_pdb_strand_id.items():
            representative_auth_chain_id = None
            if canonical_residue_list:
                 first_canon_res = canonical_residue_list[0]
                 for auth_id_3d, res3d_obj in existing_3d_residues_by_auth_id.items():
                     if res3d_obj.auth and res3d_obj.auth.number == first_canon_res['auth_seq_num'] and \
                        res3d_obj.auth.icode == first_canon_res['pdb_ins_code'] and \
                        res3d_obj.auth.chain == pdb_strand_id:
                          representative_auth_chain_id = res3d_obj.auth.chain; break
                 if not representative_auth_chain_id and pdb_strand_id in chain_molecule_types:
                     representative_auth_chain_id = pdb_strand_id
            if not representative_auth_chain_id or representative_auth_chain_id not in rna_chain_auth_ids:
                continue

            canonical_sequence = "".join([res['std_one_letter_code'] for res in canonical_residue_list])
            corrected_sequence_chars, corrected_structure_chars = [], []
            rnapolis_raw_seq, rnapolis_raw_struct = "N/A", "N/A"
            
            target_key_for_raw = representative_auth_chain_id
            if target_key_for_raw in rnapolis_parsed_data:
                rnapolis_raw_seq, rnapolis_raw_struct = rnapolis_parsed_data[target_key_for_raw]
            else:
                found_raw = False
                for rnap_key_iter in rnapolis_parsed_data.keys():
                    if rnap_key_iter.upper() == target_key_for_raw.upper():
                        rnapolis_raw_seq, rnapolis_raw_struct = rnapolis_parsed_data[rnap_key_iter]
                        found_raw = True; break
                if not found_raw and pdb_strand_id in rnapolis_parsed_data:
                     rnapolis_raw_seq, rnapolis_raw_struct = rnapolis_parsed_data[pdb_strand_id]

            for canon_res in canonical_residue_list:
                auth_id_to_lookup = (representative_auth_chain_id, canon_res['auth_seq_num'], canon_res['pdb_ins_code'])
                res3d = existing_3d_residues_by_auth_id.get(auth_id_to_lookup)
                if res3d is None:
                    corrected_sequence_chars.append('-'); corrected_structure_chars.append('.')
                else:
                    rnapolis_idx = residue_to_rnapolis_idx.get(res3d)
                    if rnapolis_idx is None:
                        corrected_sequence_chars.append(canon_res['std_one_letter_code']); corrected_structure_chars.append('.')
                    else:
                        _, struct_char = rnapolis_idx_to_chars.get(rnapolis_idx, ('?', '.'))
                        corrected_sequence_chars.append(canon_res['std_one_letter_code']); corrected_structure_chars.append(struct_char)
            
            corrected_sequence = "".join(corrected_sequence_chars)
            corrected_structure = "".join(corrected_structure_chars)
            structure_type_str = rna_structure_type.get(representative_auth_chain_id, "N/A")
            partners = pairing_partners.get(representative_auth_chain_id, set())
            other_partners = partners - {representative_auth_chain_id}
            partners_str = ", ".join(sorted(list(other_partners))) if other_partners else ("无" if structure_type_str == "intramolecular" else "N/A")
            if structure_type_str == "no_secondary_structure_found": partners_str = "无"

            file_rows.append([
                pdb_id, pdb_strand_id, canonical_sequence or "N/A",
                rnapolis_raw_seq, rnapolis_raw_struct,
                corrected_sequence or "N/A", corrected_structure or "N/A",
                resolution, structure_type_str, partners_str
            ])
        logger.info(f"完成文件: {pdb_id} / Finished file: {pdb_id}")

    except Exception as e:
        logger.error(f"处理文件 {pdb_id} 时发生主错误: {e} / "
                     f"Error occurred while processing file {pdb_id}: {e}", exc_info=False)
        return []
    finally:
        if temp_cif_handle and hasattr(temp_cif_handle, 'close') and not temp_cif_handle.closed:
             try: temp_cif_handle.close()
             except Exception as e_close:
                 logger.warning(f"无法关闭临时文件句柄 for {pdb_id}: {e_close} / "
                                f"Cannot close temp file handle for {pdb_id}: {e_close}")
    return file_rows

def run_batch_analysis(input_dir: str, output_csv: str, mod_cache_path: str, log_filepath: str, max_workers: int = None):
    """
    批量分析指定目录下的 mmCIF 文件，输出 RNA 链分析结果
    Batch analyze mmCIF files in a directory and output RNA chain analysis results.

    Args:
        input_dir (str): mmCIF 文件目录 / Directory containing mmCIF files.
        output_csv (str): 输出 CSV 路径 / Output CSV file path.
        mod_cache_path (str): modifications_cache.json 路径 / Path to modifications_cache.json.
        log_filepath (str): 日志文件路径 / Log file path.
        max_workers (int, optional): 并行进程数 / Number of parallel workers.
    """
    log_level_str = os.getenv("LOGLEVEL", "INFO")
    setup_logging(log_level_str, log_filepath)
    logging.getLogger("rnapolis").setLevel(logging.WARNING)

    modification_handler = None
    if not os.path.exists(mod_cache_path):
        logger.error(f"错误：modifications_cache.json 文件未找到于: {mod_cache_path} / "
                     f"Error: modifications_cache.json not found at: {mod_cache_path}")
        modification_handler = ModificationHandler()
    else:
        if USING_REAL_MODIFICATION_HANDLER:
            try:
                modification_handler = ModificationHandler(mod_cache_path)
                logger.info(f"成功加载 ModificationHandler 从: {mod_cache_path} / "
                            f"Successfully loaded ModificationHandler from: {mod_cache_path}")
            except Exception as e:
                logger.error(f"加载 ModificationHandler 时出错: {e} / "
                             f"Error loading ModificationHandler: {e}", exc_info=True)
                modification_handler = ModificationHandler()
        else:
            modification_handler = ModificationHandler()

    start_time = time.time()
    logger.info(f"脚本启动: RNA 结构分析 (RDL1.6.7-1.4) / Script started: RNA structure analysis (RDL1.6.7-1.4)")
    logger.info(f"输入目录: {input_dir} / Input directory: {input_dir}")
    logger.info(f"输出 CSV: {output_csv} / Output CSV: {output_csv}")
    logger.info(f"日志文件: {log_filepath} / Log file: {log_filepath}")

    try:
        all_files = [f for f in os.listdir(input_dir) if f.lower().endswith((".cif", ".cif.gz"))]
        cif_files = [os.path.join(input_dir, f) for f in all_files]
        if not cif_files:
            logger.warning(f"目录 {input_dir} 中未找到 .cif 或 .cif.gz 文件。/ "
                           f"No .cif or .cif.gz files found in directory {input_dir}.")
            return
        logger.info(f"找到 {len(cif_files)} 个文件待处理。/ Found {len(cif_files)} files to process.")
    except FileNotFoundError:
        logger.error(f"输入目录不存在: {input_dir} / Input directory does not exist: {input_dir}")
        return
    except Exception as e:
        logger.error(f"读取输入目录时出错: {e} / Error reading input directory: {e}")
        return

    try:
        with open(output_csv, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([
                'PDB_ID', 'Chain_ID', 'Sequence_Canonical',
                'Sequence_RNApolis_Raw', 'Structure_RNApolis_Raw',
                'Sequence_Corrected', 'Structure_Corrected',
                'Resolution', 'Pairing_Type', 'Pairing_Partners'
            ])

            num_workers = max_workers if max_workers else os.cpu_count()
            logger.info(f"使用 {num_workers} 个进程进行并行处理... / Using {num_workers} processes for parallel processing...")

            processed_files_count = 0
            error_files_count = 0

            with ProcessPoolExecutor(max_workers=num_workers) as executor:
                process_file_with_handler = partial(process_file, modification_handler=modification_handler)
                futures = {executor.submit(process_file_with_handler, cif_path): cif_path for cif_path in cif_files}

                for future in tqdm(as_completed(futures), total=len(cif_files), desc="处理 CIF 文件 / Processing CIF files"):
                    cif_path = futures[future]
                    pdb_id_for_log = os.path.splitext(os.path.basename(cif_path))[0]
                    try:
                        file_result_rows = future.result()
                        if file_result_rows:
                            writer.writerows(file_result_rows)
                            processed_files_count +=1
                    except Exception as exc:
                        logger.error(f"文件 {pdb_id_for_log} 在子进程中处理失败: {exc} / "
                                     f"File {pdb_id_for_log} failed in subprocess: {exc}", exc_info=False)
                        error_files_count += 1
            
            skipped_files_count = len(cif_files) - processed_files_count - error_files_count

    except IOError as e:
        logger.error(f"无法写入 CSV 文件 {output_csv}: {e} / "
                     f"Cannot write to CSV file {output_csv}: {e}")
        return
    except Exception as e:
        logger.error(f"主程序发生未知错误: {e} / Unknown error occurred in main program: {e}", exc_info=True)
        return

    end_time = time.time()
    duration = end_time - start_time
    logger.info("-" * 50)
    logger.info(f"分析完成! / Analysis complete!")
    logger.info(f"总共检查文件数: {len(cif_files)} / Total files checked: {len(cif_files)}")
    logger.info(f"成功处理并写入数据的文件数: {processed_files_count} / Files successfully processed and written: {processed_files_count}")
    logger.info(f"处理失败的文件数 (子进程错误): {error_files_count} / Files failed (subprocess errors): {error_files_count}")
    logger.info(f"跳过的文件数 (无RNA链或RNApolis无结果等): {skipped_files_count} / Files skipped (no RNA chain or no RNApolis result): {skipped_files_count}")
    logger.info(f"总耗时: {duration:.2f} 秒 / Total time: {duration:.2f} seconds")
    logger.info(f"结果已保存至 CSV: {output_csv} / Results saved to CSV: {output_csv}")
    logger.info(f"详细日志已保存至: {log_filepath} / Detailed log saved to: {log_filepath}")
    logger.info("-" * 50)

def main():
    """
    主入口，解析命令行参数并启动分析
    Main entry: parse command-line arguments and start analysis.
    """
    parser = argparse.ArgumentParser(
        description="Analyze RNA chains in mmCIF files and output a CSV report. "
                    "分析 mmCIF 文件中的 RNA 链并输出 CSV 报告。"
    )
    parser.add_argument("--input_dir", type=str, default="./data/downloaded_mmcif",
                        help="Directory containing mmCIF files. 包含 mmCIF 文件的目录。")
    parser.add_argument("--output_csv", type=str, default="./data/analysis_output/02_rna_chain_analysis.csv",
                        help="Path to output CSV file. 输出 CSV 文件路径。")
    parser.add_argument("--mod_cache", type=str, default="./utils/modifications_cache.json",
                        help="Path to modifications_cache.json. modifications_cache.json 路径。")
    parser.add_argument("--log_file", type=str, default="./data/analysis_output/02_analyze_rna_chains.log",
                        help="Path to log file. 日志文件路径。")
    parser.add_argument("--max_workers", type=int, default=80,
                        help="Number of parallel workers (default: 8). 并行进程数（默认8）。")

    args = parser.parse_args()

    run_batch_analysis(
        input_dir=args.input_dir,
        output_csv=args.output_csv,
        mod_cache_path=args.mod_cache,
        log_filepath=args.log_file,
        max_workers=args.max_workers
    )

if __name__ == "__main__":
    main()