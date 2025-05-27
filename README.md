# RNA3DB-2D-Structures: RNA3DB RNA Secondary Structure Dataset Generation Pipeline

This project provides a comprehensive pipeline to process 3D RNA structures from the Protein Data Bank (PDB), as cataloged by RNA3DB, to generate a curated dataset of RNA secondary structures in BPSEQ format. The pipeline automates downloading mmCIF files, analyzing RNA chains using RNApolis, filtering and selecting representative chains, re-splitting datasets, and finally, converting them into a 2D structure format.

The workflow is primarily based on the methodologies and scripts developed in the RDL1.6.6 and RDL1.6.7 series, adapted for robust execution and ease of use.

## Features

* **Automated Workflow**: End-to-end pipeline script (`run_pipeline.sh`) for one-click execution.
* **Data Ingestion**: Downloads mmCIF files based on PDB IDs from an RNA3DB `cluster.json`.
* **RNA Chain Analysis**: Utilizes `rnapolis-py` for detailed 3D and 2D structure analysis of RNA chains.
* **Canonicalization**: Standardizes modified RNA residues to A, C, G, U, or N.
* **Structure Correction**: Aligns RNApolis-derived secondary structures with canonical PDB sequences.
* **Filtering & Selection**: Implements logic to filter chains (e.g., intramolecular pairing, sequence coverage) and select the best representative chain per RNA3DB cluster based on resolution and coverage.
* **Dataset Splitting**: Leverages `rna3db`'s splitting functionality to create train/validation/test sets from selected representative clusters.
* **BPSEQ Output**: Generates final secondary structure data in the widely used BPSEQ format, organized by dataset.
* **Customizable**: Individual scripts can be run with different parameters, and the pipeline configuration is flexible.

## Workflow Overview

The pipeline consists of the following main steps, orchestrated by `run_pipeline.sh`:

```mermaid
graph TD
    A[Input: RNA3DB cluster.json] --> B(Step 1: Download Raw mmCIFs `01_download_raw_mmcif.py`);
    B --> C(Step 2: Analyze RNA Chains `02_analyze_rna_chains.py`);
    C --> D(Output: rna_chain_analysis.csv);
    A --> E(Step 3: Filter & Select Chains `03_filter_select_chains.py`);
    D --> E;
    E --> F(Output: filtered_selected_chains.csv);
    E --> G(Output: selected_clusters.json);
    G --> H(Step 4: Resplit Dataset `04_resplit_dataset.sh`);
    H --> I(Output: resplit_dataset.json);
    F --> J(Step 5: Assign Dataset Info `05_assign_dataset_info.py`);
    I --> J;
    J --> K(Output: final_selected_seqs_info.csv);
    K --> L(Step 6: Generate BPSEQ Files `06_generate_bpseq.py`);
    L --> M[Output BPSEQ files];  

    subgraph "User Provided Inputs"
        A
    end
    subgraph "Pipeline Outputs"
        M
    end

```

1. **Download Raw mmCIFs (`01_download_raw_mmcif.py`)**:
   * Reads PDB IDs from the user-provided `cluster.json`.
   * Uses `utils/batch_download.sh` to download corresponding mmCIF files (as `.cif.gz`) from RCSB.
   * Uncompresses the downloaded files to `.cif`.
2. **Analyze RNA Chains (`02_analyze_rna_chains.py`)**:
   * Processes each downloaded `.cif` file.
   * Uses `rnapolis-py` to identify RNA chains, parse 3D structure, and predict 2D structure (dot-bracket).
   * Extracts canonical sequences from `_pdbx_poly_seq_scheme`.
   * Generates "corrected" sequences and dot-bracket structures aligned to the canonical sequence length.
   * Outputs a detailed CSV (`rna_chain_analysis.csv`) with information for each RNA chain.
3. **Filter & Select Chains (`03_filter_select_chains.py`)**:
   * Loads `rna_chain_analysis.csv` and the original `cluster.json`.
   * Filters chains (e.g., intramolecular pairing, minimum sequence coverage from RNApolis).
   * Merges analysis data with metadata from `cluster.json`.
   * For each representative PDB cluster (`repr_pdb_chain` from `cluster.json`), selects the "best" member chain based on resolution and sequence coverage.
   * Outputs `filtered_selected_chains.csv` (all chains passing filters, with a selection mark) and `selected_clusters.json` (a subset of the original `cluster.json` containing only the representative clusters that have a selected member).
4. **Resplit Dataset (`04_resplit_dataset.sh`)**:
   * Takes `selected_clusters.json` as input.
   * Uses the `rna3db split` command-line tool to partition these selected representative clusters into `train_set`, `valid_set`, and `test_set`.
   * Outputs `resplit_dataset.json`.
5. **Assign Dataset Info (`05_assign_dataset_info.py`)**:
   * Loads `filtered_selected_chains.csv` and `resplit_dataset.json`.
   * Assigns each chain in `filtered_selected_chains.csv` to a dataset (train/valid/test) based on its `repr_pdb_chain`'s assignment in `resplit_dataset.json`.
   * Outputs `final_all_info.csv` (all filtered chains with dataset info) and `final_selected_seqs_info.csv` (only the selected representative chains with dataset info).
6. **Generate BPSEQ Files (`06_generate_bpseq.py`)**:
   * Loads `final_selected_seqs_info.csv`.
   * Converts the `Sequence_Canonical` and `Structure_Corrected` for each selected chain into BPSEQ format.
   * Saves `.bpseq` files into subdirectories named after their assigned dataset (e.g., `bpseq_files/train_set/`).
   * Optionally creates another set of BPSEQ files for sequences shorter than 600nt.

## Prerequisites

1. **Bash Environment**: Linux or macOS (Windows users can use WSL).
2. **Python 3.8+**: With `pip` installed.
3. **External Tools**:
   * `curl`: Required by `utils/batch_download.sh`.
   * `gunzip`: Required by `01_download_raw_mmcif.py` for uncompressing files.
   * `rna3db` Python package: Required by `04_resplit_dataset.sh` (for `python -m rna3db split`) and potentially by `02_analyze_rna_chains.py` if using the official `ModificationHandler`.
     * Install via pip: `pip install rna3db`
     * Or, for the latest version, install from its GitHub repository: `pip install git+https://github.com/marcellszi/rna3db.git`
4. **Python Packages**: Listed in `requirements.txt`.

## Installation

1. **Clone the repository**:

   ```bash
   git clone https://github.com/your-username/RNA3DB-2D-Structures.git
   cd RNA3DB-2D-Structures
   ```
2. **Install Python dependencies**:
   It's highly recommended to use a virtual environment:

   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install -r requirements.txt
   ```
3. **Install `rna3db`**:
   As mentioned in prerequisites, ensure `rna3db` is installed:

   ```bash
   pip install rna3db
   # or for the latest version:
   # pip install git+https://github.com/marcellszi/rna3db.git
   ```

   Verify installation by running `python -m rna3db --help`.
4. **Ensure `batch_download.sh` is executable**:

   ```bash
   chmod +x ./utils/batch_download.sh
   ```

## Data Preparation

1. **`cluster.json`**:

   * You need an RNA3DB `cluster.json` file. This file defines the RNA 3D structure clusters.
   * You can obtain this file from an official RNA3DB release. For example, from the [RNA3DB Releases page](https://github.com/marcellszi/rna3db/releases), download the `rna3db-jsons.tar.gz` archive corresponding to a release (e.g., `2024-12-04-full-release`).
   * Extract `cluster.json` from the archive.
   * Place this `cluster.json` file in a location accessible by the pipeline (e.g., the project's `data/` directory or a custom path).
   * An `example_cluster.json` is provided in the `data/` directory for structure reference, but it's likely too small for a full run.
2. **`modifications_cache.json`**:

   * This file is used by `02_analyze_rna_chains.py` to map modified RNA residues to standard ones.
   * A version of this file, sourced from the `rna3db` library's test data, is already included in this project at `data/modifications_cache.json`. The pipeline is configured to use this by default.

## Running the Pipeline

The entire pipeline can be run using the `run_pipeline.sh` script.

1. **Configure the pipeline**:

   * Copy the configuration template:
     ```bash
     cp config/pipeline_config_template.sh config/pipeline_config.sh
     ```
   * Edit `config/pipeline_config.sh` to set the correct paths, especially:
     * `CLUSTER_JSON_PATH`: Path to your `cluster.json` file.
     * Other output directories if you wish to change them from the defaults.
     * Adjust `ANALYZE_MAX_WORKERS`, `TRAIN_RATIO`, `VALID_RATIO` as needed.
2. **Execute the pipeline script**:

   ```bash
   bash ./run_pipeline.sh
   ```

   The script will:

   * Load configurations from `config/pipeline_config.sh`.
   * Create necessary output directories.
   * Execute steps 01 through 06 in sequence.
   * Log progress and errors to the console and to `pipeline.log` (by default).

   Output files will be generated in the directories specified in `config/pipeline_config.sh` (e.g., `./downloaded_mmcif/`, `./analysis_output/`, `./bpseq_files/`).

## Running Individual Scripts

You can also run individual scripts if you need to perform specific steps or debug. Each script in the `scripts/` directory accepts command-line arguments. Use the `-h` or `--help` flag for details on each script's usage.

**Example: Running `02_analyze_rna_chains.py` manually**

```bash
python ./scripts/02_analyze_rna_chains.py \
    ./downloaded_mmcif/ \
    ./analysis_output/my_custom_analysis.csv \
    --mod_cache ./data/modifications_cache.json \
    --log_file ./analysis_output/analyze_manual.log \
    --max_workers 4
```

## Output Structure

* **`downloaded_mmcif/`**: Contains the downloaded and uncompressed `.cif` files.
* **`analysis_output/`**:
  * `rna_chain_analysis.csv`: Detailed analysis results from `02_analyze_rna_chains.py`.
  * `filtered_selected_chains.csv`: Output from `03_filter_select_chains.py` with all chains passing initial filters and selection marks.
  * `selected_clusters.json`: Subset of `cluster.json` used for re-splitting, from `03_filter_select_chains.py`.
  * `resplit_dataset.json`: Dataset split information from `04_resplit_dataset.sh`.
  * `final_all_info.csv`: All filtered chains with their assigned dataset.
  * `final_selected_seqs_info.csv`: Only the selected representative chains with their assigned dataset (input for BPSEQ generation).
  * Log files for various steps.
* **`bpseq_files/`**:
  * `train_set/`: BPSEQ files for the training set.
  * `valid_set/`: BPSEQ files for the validation set.
  * `test_set/`: BPSEQ files for the test set.
  * (Optional) `bpseq_less600nt/`: If enabled, contains BPSEQ files for sequences < 600nt, also organized by dataset.
* **`pipeline.log`**: Main log file for `run_pipeline.sh`.

## Notes on RNA3DB Core Dataset Generation

This project focuses on processing an existing RNA3DB `cluster.json` to generate a 2D structure dataset. If you need to:

* Generate the `cluster.json` file from scratch (i.e., parse all of PDB, filter, perform sequence and structural clustering as defined by the RNA3DB methodology).
* Understand the detailed methodology behind RNA3DB's clustering and component definition.

Please refer to the official **[RNA3DB GitHub repository](https://github.com/marcellszi/rna3db)** and its documentation. The `rna3db` Python package provides the necessary tools for these tasks.

## Contributing

Contributions, bug reports, and feature requests are welcome! Please open an issue or submit a pull request.

## License

This project is licensed under the [MIT License](LICENSE).

```

```
