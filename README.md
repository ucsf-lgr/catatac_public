># Project README

CATATAC Public Repo

## Directory Structure

### dasatinib_primary_01HD/
Contains primary analysis data organized by experimental conditions:
- **DMSO1/** - Control condition replicate 1
- **DMSO2/** - Control condition replicate 2
- **DASA1/** - Dasatinib treatment replicate 1
- **DASA2/** - Dasatinib treatment replicate 2

**Inputs**: FASTQS, library files, protospacer information files, path to CATATAC primary pipeline
**Outputs**: Cellranger output, csvs

### secondary/
Contains secondary analysis notebooks, organized by analysis step numbers.

#### Data Processing and QC Notebooks
- Secondary CATATAC preprocessing pipeline
	- **DMSO1/** - Control condition replicate 1
	- **DMSO2/** - Control condition replicate 2
	- **DASA1/** - Dasatinib treatment replicate 1
	- **DASA2/** - Dasatinib treatment replicate 2
	- **Inputs**: Outputs from primary, path to CATATAC secondary pipeline
	- **Outputs**: RDS, sobj, csvs, guide caller calls, plots

- **06_filter_cc_genes_and_mixscape_all_samples_ngc.ipynb**
  - **Description**: Filters cell cycle genes and performs Mixscape analysis. Assigns calls
  - **Inputs**: sobj/RDS, guide caller calls
  - **Outputs**: sobj/RDS, plots

- **07_merge_samples_ngc.ipynb**
  - **Description**: Merges multiple sample datasets for analysis.
  - **Inputs**: sobj/RDS
  - **Outputs**: sobj/RDS, csv, plots

#### Motif Analysis
- **07b_ngc_motif_analysis.ipynb**
  - **Description**: Performs motif analysis.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

- **07c_ngc_motif_footprinting.ipynb**
  - **Description**: Analyzes motif footprinting.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

#### Gene Regulatory Network Analysis
- **08_output_notebook.ipynb**
  - **Description**: Contains output results from previous analyses.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

- **08b_celloracle.ipynb**
  - **Description**: Applies CellOracle for gene regulatory network analysis.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

- **08ba_cicero.R.ipynb**
  - **Description**: Uses Cicero for cis-regulatory interaction analysis.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

- **08bb_cicero_downstream_analysis.R.ipynb**
  - **Description**: Performs downstream analysis of Cicero results.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

- **08bc_cicero_link_plots_for_ks.R.ipynb**
  - **Description**: Generates plots for Cicero links.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

- **08c_pando.ipynb**, **08c_pando_ngc.ipynb**, **08c_pando_subsets.ipynb**
  - **Description**: Applies Pando for regulatory network inference.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

#### Visualization and Figure Generation
- **15_pseudobulk_atac_heatmaps.R.ipynb**
  - **Description**: Creates heatmaps from pseudobulk ATAC-seq data.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

- **16_promoter_expression_plots.R.ipynb**
  - **Description**: Generates plots of promoter expression.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

- **17_make_figure2b.ipynb**, **17_make_figure2b_output.ipynb**, **17b_figure2b_remake.R.ipynb**
  - **Description**: Creates Figure 2B.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

- **18_zfmp2_coef_heatmaps.R.ipynb**
  - **Description**: Generates heatmaps of ZFPM2 coefficients.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

- **20_make_aggregate_heatmap.R.ipynb**
  - **Description**: Creates aggregate heatmaps.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

- **22_make_new_dotplots.R.ipynb**
  - **Description**: Generates dotplot visualizations.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

#### Specific Analyses
- **19_check_gata2_zfpm2_overlap.R.ipynb**
  - **Description**: Analyzes overlap between GATA2 and ZFPM2 regions.
  - **Inputs**: CSV file with coefficients containing TF binding data
  - **Outputs**: CSV file with overlapping regions between GATA2 and ZFPM2 ("19_overlapping_regions_gata2_zfpm2.csv")

- **21_all_zfpm2_coefs_atac_heatmap.R.ipynb**
  - **Description**: Creates heatmap of ZFPM2 coefficients from ATAC-seq data.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

- **23_calc_guide_knockdown.R.ipynb**
  - **Description**: Calculates guide RNA knockdown efficiency.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

- **protospacer_distribution_plot.ipynb**
  - **Description**: Plots distribution of protospacers.
  - **Inputs**: _[To be filled in by user]_
  - **Outputs**: _[To be filled in by user]_

### secondary/fimo/
- **Description**: _[To be filled in by user]_
- **Inputs**: _[To be filled in by user]_
- **Outputs**: _[To be filled in by user]_

### secondary/gata1_filter/
- **Description**: _[To be filled in by user]_
- **Inputs**: _[To be filled in by user]_
- **Outputs**: _[To be filled in by user]_

### secondary/more_atac_analysis/
- **Description**: _[To be filled in by user]_
- **Inputs**: _[To be filled in by user]_
- **Outputs**: _[To be filled in by user]_

### secondary/23_make_qc_comparisons/
- **Description**: _[To be filled in by user]_
- **Inputs**: _[To be filled in by user]_
- **Outputs**: _[To be filled in by user]_

## Workflow Summary

The overall workflow appears to follow these general steps:

1. **Primary Analysis**: Processing raw data in the dasatinib_primary_01HD directory, organized by experimental conditions.
2. **Data Processing and QC**: Filtering, quality control, and preprocessing of the data (06_, 07_ notebooks).
3. **Motif Analysis**: Analysis of DNA binding motifs (07b_, 07c_ notebooks).
4. **Gene Regulatory Network Analysis**: Inference of gene regulatory networks using various methods like CellOracle, Cicero, and Pando (08_ notebooks).
5. **Specific Analyses**: Focused analyses on particular transcription factors such as GATA2 and ZFPM2 (19_, 21_ notebooks).
6. **Visualization**: Generation of figures and plots for publication or interpretation (15_, 16_, 17_, 18_, 20_, 22_ notebooks).

## Common File Inputs/Outputs

Most notebooks in this project typically:

- **Read**: 
  - Seurat objects (.rds files)
  - Count matrices (CSV, TSV, or HDF5 format)
  - Peak files (BED format)
  - Annotation files (GTF/GFF)
  - Previously generated results (CSV, RDS)

- **Write**:
  - Processed data (CSV, RDS)
  - Analysis results (CSV, TSV)
  - Plots and figures (PDF, PNG)

---