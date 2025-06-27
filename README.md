# AFChimera - AlphaFold for chimeric proteins with windowed MSA

Implements windowed MSA that concatenates multiple sequence alignments (MSAs) for structure prediction of chimeric proteins.

## Installation

AFChimera requires only Python 3.6+ with standard library modules. No additional dependencies needed.

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

## Citation

If you use AFChimera in your research, please cite:

```
@article{afchimera2024,
  title={Improving Prediction Accuracy in Chimeric Proteins with Windowed Multiple Sequence Alignment},
  author={Sanketh Vedula and Alex M. Bronstein and Ailie Marx},
  year={2025},
}

```