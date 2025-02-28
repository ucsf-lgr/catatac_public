# CATATAC Public Repo

## Directory Structure

### dasatinib_primary_01HD/
Contains primary analysis data organized by experimental conditions:
- **DMSO1/** - Control condition replicate 1
- **DMSO2/** - Control condition replicate 2
- **DASA1/** - Dasatinib treatment replicate 1
- **DASA2/** - Dasatinib treatment replicate 2

- **Use**: Seperate pipeline README in 00_catatac_pipeline
- **Inputs**: FASTQS, library files, protospacer information files, path to CATATAC primary pipeline
- **Outputs**: Cellranger output, csvs, tsvs

### secondary/
Contains secondary analysis notebooks, organized by analysis step numbers.

#### Data Processing and QC Notebooks
- Secondary CATATAC preprocessing pipeline
	- **DMSO1/** - Control condition replicate 1
	- **DMSO2/** - Control condition replicate 2
	- **DASA1/** - Dasatinib treatment replicate 1
	- **DASA2/** - Dasatinib treatment replicate 2
  - **Use**: Seperate pipeline README in 00_catatac_pipeline
	- **Inputs**: Outputs from primary, path to CATATAC secondary pipeline
	- **Outputs**: RDS, sobj, csvs, guide caller calls, plots, notebooks for each step in the pipeline

- **00_filter_cc_genes_and_mixscape_all_samples_supp3C.ipynb**
  - **Description**: Filters cell cycle genes and performs Mixscape analysis. Assigns calls
  - **Inputs**: sobj/RDS, guide caller calls
  - **Outputs**: sobj/RDS, plots, Supplemental Figure 3C

- **01_merge_samples_fig3EF.ipynb**
  - **Description**: Merges multiple sample datasets for analysis.
  - **Inputs**: sobj/RDS
  - **Outputs**: sobj/RDS, csv, plots, Figures 3 E and F

#### Motif Analysis
- **02_ngc_motif_analysis.ipynb**
  - **Description**: Performs motif analysis.
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: Motif csvs

#### Gene Regulatory Network Analysis
- **03_pando_ngc.R**
  - **Description**: Applies Pando for regulatory network inference.
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: GRN sobj and modules

- **04_pando_ngc_fig4D_supp4.ipynb**
  - **Description**: Makes plots with outputs from 03.
  - **Inputs**: Merged sobj/RDS + 03 outs
  - **Outputs**: Figure 4D and Supplemental Figure 4, GRN coefs

#### Visualization and Figure Generation 
- **05_DEA_heatmaps_fig4A.R.ipynb**
  - **Description**: Creates heatmaps from DEA data.
  - **Inputs**: Merged sobj/RDS, DEA csvs
  - **Outputs**: Plot Fig 4A

- **06_promoter_expression_plots_supp3B.R.ipynb**
  - **Description**: Generates plots of promoter expression. 
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: Supplemental Fig 3B plot

- **07_fig4EF.R.ipynb**
  - **Description**: Creates Fig 4E and F.
  - **Inputs**: Merged sobj/RDS, GRN coefs (saved in 04)
  - **Outputs**: Figures 4E and F

- **08_fig3C.R.ipynb**
  - **Description**: Generates dotplot visualizations showing target gene knockdown in cell populations which received the indicated guide in both dasatinib and DMSO treated conditions.
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: Figure 3C

- **09_qc_multiomic_comparisons_supp1CDEF.R**
  - **Description**: Generates QC comparison plots across methods. [Link to paper where external data is from](https://www.sciencedirect.com/science/article/pii/S2405471224003661?via%3Dihub)
  - **Inputs**: Merged sobj/RDS, GSM8528725_MPS, CROP-Multiome.RDS, qc_Source_Data
  - **Outputs**: Supplemental Figures 1C , D, E, and F

#### Specific Analyses
- **10_calc_guide_knockdown.R.ipynb**
  - **Description**: Calculates guide RNA knockdown efficiency.
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: Table of knockdown

- **11_protospacer_distribution_plot_supp3A.ipynb**
  - **Description**: Plots distribution of protospacers.
  - **Inputs**: Guide calls
  - **Outputs**: Supplemental Figure 3A

## Workflow Summary

The overall workflow appears to follow these general steps:

1. **Primary Analysis**: Processing raw data in the dasatinib_primary_01HD directory, organized by experimental conditions.
2. **Data Processing and QC**: Filtering, quality control, and preprocessing of the data (06_, 07_ notebooks).
3. **Motif Analysis**: Analysis of DNA binding motifs (07b_, 07c_ notebooks).
4. **Gene Regulatory Network Analysis**: Inference of gene regulatory networks using Pando (08_ notebooks).
5. **Specific Analyses**: Focused on guide knockdown and protospacer distribution.
6. **Visualization**: Generation of figures and plots for publication or interpretation.

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