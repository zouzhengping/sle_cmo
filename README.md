# TSCM Signature Analysis in SLE - Analysis Code

This repository contains the analysis code for the manuscript: **"Identification of Stem-like Memory T Cells in Systemic Lupus Erythematosus and Their Association with Disease Activity"**

## Code Structure

### Main Analysis Pipeline

The analysis follows a staged pipeline approach:

#### Stage 4: Single-cell RNA-seq Analysis (GSE137029)
- `40_stage4_basic_qc_sle.py` - Basic quality control
- `41_stage4_filter_cells_sle.py` - Cell filtering
- `42_stage4_normalization_and_hvg.py` - Normalization and highly variable gene identification
- `43_stage4_tcell_score_and_subset_sle.py` - T-cell identification and subsetting
- `44_stage4_tcells_cluster_only_sle.py` - T-cell clustering
- `45_stage4_tscm_scoring_sle.py` - TSCM scoring
- `46_stage4_mark_tscm_clusters.py` - TSCM cluster identification
- `47_stage4_de_tscm_vs_other.py` - Differential expression analysis
- `48_stage4_visualization.py` - Visualization of results

#### Stage 5: Signature Analysis and Enrichment
- `50_stage5A_build_pseudobulk_tscm.py` - Pseudobulk construction
- `51_stage5B_enrichment_tscm.py` - Pathway enrichment analysis (GO, KEGG, Reactome)
- `52_stage5C_ligand_receptor_tscm.py` - Ligand-receptor analysis
- `53_stage5D_trajectory_tscm.py` - Trajectory analysis

### Validation and External Datasets

#### GSE88884 Validation
- `download_gse88884_data.py` - Download GSE88884 data
- `process_gse88884_data.py` - Process GSE88884 data
- `analyze_tscm_vs_sledai.py` - TSCM signature vs SLEDAI correlation analysis
- `compare_tscm_by_sledai_groups.py` - Group comparison by SLEDAI

#### Other Validation Datasets
- `download_gse211700_data.py` - Download GSE211700 data
- `validate_gse211700_tscm_signature.py` - Validate signature in GSE211700
- `download_gse61635_data.py` - Download GSE61635 data
- `validate_gse61635_tscm_signature.py` - Validate signature in GSE61635
- `convert_gse61635_probe_to_gene.R` - R script for probe-to-gene conversion

### Signature Refinement
- `create_core_tscm_signature.py` - Create core TSCM signature (30-50 genes)
- `create_signature_comparison.py` - Compare core vs full signature

### Figure Generation
- `create_schematic_figure.py` - Generate workflow schematic (Figure S1)
- `fix_and_regenerate_figures.py` - Fix and regenerate figures for PLOS ONE format
- `optimize_figures_for_plos_one.py` - Optimize figures for publication
- `convert_figures_for_plos_one.py` - Convert figures to PLOS ONE format

### Utility Scripts
- `deconvolution_tscm_proportion.py` - Deconvolution analysis
- `validate_with_tcell_control.py` - Validation with T-cell controls
- `collect_references.py` - Collect references
- `check_sledai_metadata.py` - Check SLEDAI metadata

## Dependencies

### Python Packages
- scanpy >= 1.9.0
- pandas >= 1.5.0
- numpy >= 1.23.0
- scipy >= 1.9.0
- matplotlib >= 3.6.0
- seaborn >= 0.12.0
- scikit-learn >= 1.1.0
- statsmodels >= 0.13.0
- gseapy >= 1.0.0
- anndata >= 0.8.0

### R Packages (for some scripts)
- GEOquery
- limma
- Biobase

## Usage

### Running the Main Pipeline

1. **Download data**: First download the required datasets (GSE137029, GSE88884, etc.)
2. **Run Stage 4 scripts in order**: Execute scripts 40-48 sequentially
3. **Run Stage 5 scripts**: Execute scripts 50-53 for enrichment analysis
4. **Run validation**: Execute validation scripts for external datasets

### Example

```bash
# Stage 4: Single-cell analysis
python scripts/40_stage4_basic_qc_sle.py
python scripts/41_stage4_filter_cells_sle.py
python scripts/42_stage4_normalization_and_hvg.py
# ... continue with remaining scripts

# Stage 5: Enrichment analysis
python scripts/50_stage5A_build_pseudobulk_tscm.py
python scripts/51_stage5B_enrichment_tscm.py
```

## Data Availability

The raw data used in this analysis are available from:
- GSE137029: Single-cell RNA-seq of SLE PBMC
- GSE88884: Bulk RNA-seq validation cohort
- GSE211700: Additional validation dataset
- GSE61635: Additional validation dataset

## Output

Results are saved in:
- `results/tables/` - Tables (TSV format)
- `results/figures/` - Figures (TIFF format for publication)
- `data/processed/` - Processed data files (H5AD format)

## Citation

If you use this code, please cite:
[Manuscript citation will be added upon publication]

## License

[Specify license]

## Contact

For questions about the code, please contact the corresponding author.

