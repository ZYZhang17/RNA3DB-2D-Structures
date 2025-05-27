"""
modification_utils.py

The ModificationHandler class in this file is inspired by and adapted from
the ModificationHandler implementation in the rna3db project (https://github.com/marcellszi/rna3db),
specifically from rna3db/parser.py.
For original implementation details, please refer to the above open-source project.

Reference:
- https://github.com/marcellszi/rna3db/blob/main/rna3db/parser.py
"""

import json
from pathlib import Path
from typing import Sequence, Union # Added Union for PathLike
import os # Added os for PathLike

# Define PathLike if not already defined elsewhere in your project's utils
PathLike = Union[str, os.PathLike]

class ModificationHandler:
    """
    Handles conversion of 3-letter residue codes (including modifications)
    to 1-letter codes, based on a pre-generated cache file.
    处理三字母残基代码（包括修饰）到单字母代码的转换，
    基于预生成的缓存文件。
    """
    # Default path relative to this file's location if no path is provided.
    # Assumes modifications_cache.json is in the same directory as this script.
    DEFAULT_CACHE_FILENAME = "modifications_cache.json"

    def __init__(self, json_path: PathLike = None):
        """
        Initializes the handler by loading the modifications cache.
        通过加载修饰缓存来初始化处理器。

        Args:
            json_path (PathLike, optional): 
                Path to `modifications_cache.json`. 
                If None, attempts to load from a default location relative to this script
                or from a few common project locations.
                `modifications_cache.json` 的路径。
                如果为 None，则尝试从相对于此脚本的默认位置或几个常见的项目位置加载。
        """
        if json_path:
            self.cache_file_path = Path(json_path)
        else:
            # Default search strategy
            # 1. Same directory as this script
            script_dir = Path(__file__).parent
            default_path_in_utils = script_dir / self.DEFAULT_CACHE_FILENAME

            # # 2. Project root/data/ (common alternative)
            # project_root_data_path = script_dir.parent / "data" / self.DEFAULT_CACHE_FILENAME

            if default_path_in_utils.is_file():
                self.cache_file_path = default_path_in_utils
            # elif project_root_data_path.is_file(): # Check data/ as a fallback if utils/ fails
            #     self.cache_file_path = project_root_data_path
                # print(f"Info: Loaded modifications_cache.json from {self.cache_file_path} (fallback to data/).")
            else:
                raise FileNotFoundError(
                    f"Could not find '{self.DEFAULT_CACHE_FILENAME}' in default locations: "
                    # f"{default_path_in_utils} or {project_root_data_path}. "
                     f"{default_path_in_utils}. "
                    f"Please provide a valid path to ModificationHandler."
                )
        
        if not self.cache_file_path.is_file():
            raise FileNotFoundError(
                f"Specified modifications_cache.json not found at: {self.cache_file_path}"
            )

        try:
            with open(self.cache_file_path, 'r') as f:
                self.modifications = json.load(f)
            # print(f"Successfully loaded modifications cache from: {self.cache_file_path}")
        except json.JSONDecodeError as e:
            raise ValueError(f"Error decoding JSON from {self.cache_file_path}: {e}")
        except Exception as e:
            raise RuntimeError(f"Could not load or parse {self.cache_file_path}: {e}")

    def is_rna(self, three_letter_code: str) -> bool:
        """
        Checks if the given 3-letter code corresponds to an RNA nucleic acid
        based on the loaded cache.
        根据加载的缓存检查给定的三字母代码是否对应于 RNA 核酸。

        Args:
            three_letter_code (str): The 3-letter residue code.
                                     三字母残基代码。
        Returns:
            bool: True if it's an RNA code, False otherwise.
                  如果是 RNA 代码则为 True，否则为 False。
        """
        return three_letter_code.upper() in self.modifications.get("rna", {})

    def is_protein(self, three_letter_code: str) -> bool:
        """
        Checks if the given 3-letter code corresponds to an amino acid
        based on the loaded cache.
        根据加载的缓存检查给定的三字母代码是否对应于氨基酸。

        Args:
            three_letter_code (str): The 3-letter residue code.
                                     三字母残基代码。
        Returns:
            bool: True if it's a protein code, False otherwise.
                  如果是蛋白质代码则为 True，否则为 False。
        """
        return three_letter_code.upper() in self.modifications.get("protein", {})

    def rna_letters_3to1(self, three_letter_code: str) -> str:
        """
        Converts a 3-letter RNA nucleic acid code to its 1-letter equivalent.
        Returns 'N' if the code is not found or not a standard RNA base.
        将三字母 RNA 核酸代码转换为其单字母等效代码。
        如果代码未找到或不是标准 RNA 碱基，则返回 'N'。

        Args:
            three_letter_code (str): The 3-letter residue code.
                                     三字母残基代码。
        Returns:
            str: The 1-letter RNA code ('A', 'C', 'G', 'U', or 'N').
                 单字母 RNA 代码 ('A', 'C', 'G', 'U', 或 'N')。
        """
        return self.modifications.get("rna", {}).get(three_letter_code.upper(), 'N')

    def protein_letters_3to1(self, three_letter_code: str) -> str:
        """
        Converts a 3-letter amino acid code to its 1-letter equivalent.
        Returns 'X' if the code is not found.
        将三字母氨基酸代码转换为其单字母等效代码。
        如果代码未找到，则返回 'X'。

        Args:
            three_letter_code (str): The 3-letter residue code.
                                     三字母残基代码。
        Returns:
            str: The 1-letter amino acid code (or 'X').
                 单字母氨基酸代码 (或 'X')。
        """
        return self.modifications.get("protein", {}).get(three_letter_code.upper(), 'X')