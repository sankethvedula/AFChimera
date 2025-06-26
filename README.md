# AFChimera - AlphaFold for chimeric proteins with windowed MSA

Implements windowed MSA that concatenates multiple sequence alignments (MSAs) in `a3m` format to prepare MSAs for structure prediction of chimeric proteins.

## Overview

AFChimera is a command-line tool that concatenates MSA files for chimeric protein structure prediction using AlphaFold. It creates windowed MSAs with N-terminal and C-terminal fusions (tag-scaffold and scaffold-tag).

The tool is optimized for preparing MSAs for AlphaFold protein structure prediction workflows.

## Installation

AFChimera requires only Python 3.6+ with standard library modules. No additional dependencies needed.

```bash
# Clone or download the AFChimera tool
# No installation required - run directly
```

## Usage

### General Usage

```bash
# Standard windowed concatenation (creates both N-terminal and C-terminal variants)
python afchimera.py --scaffold-msa scaffold.a3m --tag-msa tag.a3m

# With custom output options
python afchimera.py --scaffold-msa scaffold.a3m --tag-msa tag.a3m --out-msas-folder ./results --output-file my_chimera
```

### Parameters

- `--scaffold-msa` (required): Path to the scaffold sequence MSA file (.a3m)
- `--tag-msa` (required): Path to the tag MSA file (.a3m)
- `--out-msas-folder`: Output directory (default: ./out_msas)
- `--output-file`: Base name for output files (auto-generated if not specified)

## Output

AFChimera creates the windowed MSA, written into a3m files.

## Examples

### Example 1: Standard Windowed MSA for GFP-1ERP fusion

```bash
python afchimera.py --scaffold-msa examples/chimeric_scaffold.a3m --tag-msa examples/tag.a3m
```

This creates:
- `out_msas/N_chimeric_scaffold_tag.a3m` (N-terminal fusion: tag-scaffold)
- `out_msas/C_chimeric_scaffold_tag.a3m` (C-terminal fusion: scaffold-tag)

### Example 2: Custom output

```bash
python afchimera.py --scaffold-msa gfp.a3m --tag-msa his_tag.a3m --out-msas-folder ./results --output-file gfp_his
```

This creates files in `./results/` with names like `N_gfp_his.a3m` and `C_gfp_his.a3m`.

## License

AFChimera is provided as-is for research purposes. If this tool is used, please consider citing appropriately if used in publications.

## Citation

If you use AFChimera in your research, please cite:

```
@article{afchimera2024,
  title={AFChimera: A tool for generating chimeric protein MSAs for AlphaFold prediction},
  author={[Author1 LastName1] and [Author2 LastName2] and [Author3 LastName3]},
  year={2024},
  journal={bioRxiv},
  doi={[bioRxiv DOI]},
  url={[bioRxiv URL]}
}

```