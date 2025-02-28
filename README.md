# CATATAC Public Repo

## Directory Structure

### dasatinib_primary_01HD/
Contains primary analysis data organized by experimental conditions:
- **DMSO1/** - Control condition replicate 1
- **DMSO2/** - Control condition replicate 2
- **DASA1/** - Dasatinib treatment replicate 1
- **DASA2/** - Dasatinib treatment replicate 2

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
	- **Inputs**: Outputs from primary, path to CATATAC secondary pipeline
	- **Outputs**: RDS, sobj, csvs, guide caller calls, plots, notebooks for each step in the pipeline

- **06_filter_cc_genes_and_mixscape_all_samples_ngc.ipynb**
  - **Description**: Filters cell cycle genes and performs Mixscape analysis. Assigns calls
  - **Inputs**: sobj/RDS, guide caller calls
  - **Outputs**: sobj/RDS, plots

- **07_merge_samples_ngc.ipynb**
  - **Description**: Merges multiple sample datasets for analysis.
  - **Inputs**: sobj/RDS
  - **Outputs**: sobj/RDS, csv, plots, umaps fig 3

#### Motif Analysis
- **07b_ngc_motif_analysis.ipynb**
  - **Description**: Performs motif analysis.
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: Motif csvs

- **07c_ngc_motif_footprinting.ipynb**
  - **Description**: Analyzes motif footprinting.
  - **Inputs**: Merged sobj/RDS, DEA results 
  - **Outputs**: Footprinting plots

#### Gene Regulatory Network Analysis
- **08c_pando_ngc.R**, **08c_pando_ngc.ipynb**, **08c_pando_subsets.R**, **08c_pando_subsets.ipynb**
  - **Description**: Applies Pando for regulatory network inference. R script is the driver, ipynb for plots.
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: GRN sobj, GRN plots, Fig 4D and Supplemental Fig 4

#### Visualization and Figure Generation 
- **13_DEA_heatmaps.R.ipynb**
  - **Description**: Creates heatmaps from DEA data.
  - **Inputs**: Merged sobj/RDS, DEA csvs
  - **Outputs**: Plot Fig 4A

- **16_promoter_expression_plots.R.ipynb**
  - **Description**: Generates plots of promoter expression. 
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: Supplemental Fig 3B plot

- **17_make_figure4ef.ipynb**, **17_make_figure4ef_output.ipynb**, **17b_figure4ef_remake.R.ipynb**
  - **Description**: Creates Fig 4E and F. Different types. Remake used in paper for F.
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: Plots Fig 4E and F

- **22_make_new_dotplots.R.ipynb**
  - **Description**: Generates dotplot visualizations.
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: Plot Fig 3C

- **24_qc_multiomic_comparisons.R**
  - **Description**: Generates QC comparison plots across methods. [Link to paper where external data is from](https://www.sciencedirect.com/science/article/pii/S2405471224003661?via%3Dihub)
  - **Inputs**: Merged sobj/RDS, GSM8528725_MPS, CROP-Multiome.RDS, qc_Source_Data
  - **Outputs**: Plot Fig 3C

#### Specific Analyses
- **23_calc_guide_knockdown.R.ipynb**
  - **Description**: Calculates guide RNA knockdown efficiency.
  - **Inputs**: Merged sobj/RDS
  - **Outputs**: Table of knockdown

- **protospacer_distribution_plot.ipynb**
  - **Description**: Plots distribution of protospacers.
  - **Inputs**: Guide calls
  - **Outputs**: Supplemental Fig 3A

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