# CATATAC Public Repo

## Directory Structure

### 00_catatac_pipeline/
Contains pipeline to run 01 and 02.

### 01_catatac_pipeline_primary_templates/
Contains primary analysis scripts organized by experimental conditions:
- **DMSO1/** - Control condition replicate 1
- **DMSO2/** - Control condition replicate 2
- **DASA1/** - Dasatinib treatment replicate 1
- **DASA2/** - Dasatinib treatment replicate 2

- **Use**: Seperate pipeline README in 00_catatac_pipeline
- **Inputs**: FASTQS, library files, protospacer information files, path to CATATAC primary pipeline
- **Outputs**: Cellranger output, csvs, tsvs

### 02_catatac_pipeline_secondary_outs/
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

### 03_additional_analysis_and_figures/
#### Figure 1 Files
- **00_catatac_secondary_pipeline_outs/**
  - **Description**: Contains output notebooks from CATATAC Secondary Pipeline for fig 1 related data
  - **Inputs**: CATATAC Primary Pipeline outs. Can be ran with any template above with modified file inputs and names.
  - **Outputs**: sobj/RDS, plots, 05_demux_guides_00_demux_by_guide_fig1E.ipynb creates Figure 1E

- **01_mixcape_bothPS_fig1FG_supp1A.ipynb**
  - **Description**: Runs Mixscape.
  - **Inputs**: sobj/RDS
  - **Outputs**: sobj/RDS, csv, plots, Figures 1 F, G, and Supplemental Figure 1A

- **02_guide_rank_plot_cond9_10_fig1H.ipynb**
  - **Description**: Makes guide rank plot.
  - **Inputs**: sobj/RDS, guide calls tsvs/csvs
  - **Outputs**: Figure 1H

- **03_signac_fig1IJ.ipynb**
  - **Description**: Signac analysis for figure 1 dataset.
  - **Inputs**: sobj/RDS, guide calls tsvs/csvs
  - **Outputs**: Figure 1 I, J

#### Data Processing and QC Notebooks
- **00_filter_cc_genes_and_mixscape_all_samples_supp3C.ipynb**
  - **Description**: Filters cell cycle genes and performs Mixscape analysis. Assigns calls
  - **Inputs**: sobj/RDS, guide caller calls
  - **Outputs**: sobj/RDS, plots, Supplemental Figure 3C

- **01_merge_samples_fig3BCEF.ipynb**
  - **Description**: Merges multiple sample datasets for analysis.
  - **Inputs**: sobj/RDS
  - **Outputs**: sobj/RDS, csv, plots, Figures 3 B, C, E, and F

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

- **08_qc_multiomic_comparisons_supp1CDEF.R**
  - **Description**: Generates QC comparison plots across methods. [Link to paper where external data is from](https://www.sciencedirect.com/science/article/pii/S2405471224003661?via%3Dihub)
  - **Inputs**: Merged sobj/RDS, GSM8528725_MPS, CROP-Multiome.RDS, qc_Source_Data
  - **Outputs**: Supplemental Figures 1C , D, E, and F

#### Specific Analyses
- **9_calc_guide_knockdown.R.ipynb**
  - **Description**: Calculates guide RNA knockdown efficiency.
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: Table of knockdown

- **10_protospacer_distribution_plot_fig1D_supp3A.ipynb**
  - **Description**: Plots distribution of protospacers.
  - **Inputs**: Guide calls
  - **Outputs**: Figure 1D, Supplemental Figure 3A

- **11_Pseudobulk_DEG_heatmap_supp3DE.ipynb**
  - **Description**: Plots heatmap of pseudobulked DEG.
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: Supplemental Figures 3D and E

- **12_ATAC_analysis_fig3GHI_supp5.ipynb**
  - **Description**: Additional ATAC analysis.
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: DAR.xlsx (differentially accessible regions), Figures 3G, H, I, and Supplemental Figure 5

## Workflow Summary

The overall workflow follows these general steps:

1. **Pipeline Setup**: Setup and configuration in the 00_catatac_pipeline directory that orchestrates the primary and secondary analyses.

2. **Primary Analysis**: Processing raw FASTQ data using the pipeline templates in 01_catatac_pipeline_primary_templates, organized by experimental conditions (DMSO1, DMSO2, DASA1, DASA2).

3. **Secondary Analysis**: Further processing and quality control in 02_catatac_pipeline_secondary_outs, producing RDS/Seurat objects and guide caller outputs.

4. **Additional Analyses and Figure Generation**: In 03_additional_analysis_and_figures:
   - Data preprocessing (filtering cell cycle genes, merging samples)
   - Motif analysis for transcription factor binding sites
   - Gene regulatory network inference using Pando
   - Differential expression analysis and visualization
   - Guide RNA knockdown efficiency analysis
   - Generation of publication-quality figures and plots

This modular workflow allows for systematic analysis of CATATAC data from raw sequencing files through to biological insights and publication-ready visualizations.

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