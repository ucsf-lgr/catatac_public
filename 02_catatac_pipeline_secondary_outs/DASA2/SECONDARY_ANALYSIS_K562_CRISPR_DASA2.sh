#!/bin/bash
set -euo pipefail

txtred=$'\e[0;31m' # Red
txtgrn=$'\e[0;32m' # Green
txtylw=$'\e[0;33m' # Yellow
txtwht=$'\e[0;37m'

R --version

CONDITION="DASA2" ##### CHANGE #####
FQ_PREFIX="CRISPR-DASA2" ##### CHANGE #####
CONDITION="$FQ_PREFIX"
PRJ_NAME="CATATAC_""$CONDITION"
H5_FILE_NAME="$PRJ_NAME"".h5seurat"
H5_SINGLETS_FILE_NAME="$PRJ_NAME""_ONLY_SINGLETS.h5seurat"
SECONDARY_A_SCRIPT_PATH="path_to_secondary_scripts"
#SECONDARY_A_SCRIPT_PATH="/home/jupyter-ayang/notebooks/cat-atac_pipeline/K562/secondary/"
PRIMARY_A_OUT_PATH="path_to_primary_out_dir" ######## CHANGE ##########
SECONDARY_A_OUT_PATH="path_to_secondary_out_dir" ######## CHANGE ##########
MAIN_DATA_PATH="path_to_fastqs"
PROTOSPACER_PATH="path_to_ps_dasatinib_screen2.tsv"

PROTOSPACER_PATH_05="path_to_ps_dasatinib_screen4.tsv"

# mkdir -p "$SECONDARY_A_OUT_PATH/seurat/"
# mkdir -p "$SECONDARY_A_OUT_PATH/notebooks/"


# # STEP 1 ----------------------------------------------------------------------------
# IPYNB_FILE=01_assign_guides_to_cells.py.ipynb
# FIRST_STEP="$SECONDARY_A_SCRIPT_PATH/$IPYNB_FILE"
# OUT_PATH="$SECONDARY_A_OUT_PATH/notebooks/$IPYNB_FILE"

# echo -e "$txtgrn \n*** STEP 1: $IPYNB_FILE $txtwht"
# papermill \
#     -p library     "$FQ_PREFIX" \
#     -p fastq_path  "$MAIN_DATA_PATH" \
#     -p tsv_path    "$PRIMARY_A_OUT_PATH/tsv" \
#     "$FIRST_STEP" \
#     "$OUT_PATH"


# # STEP 2 -----------------------------------------------------------------------------
# ##### CHANGE DATA PATH #####
# IPYNB_FILE="02_seurat.R.ipynb"
# SECOND_STEP="$SECONDARY_A_SCRIPT_PATH"/"$IPYNB_FILE"
# OUT_PATH="$SECONDARY_A_OUT_PATH/notebooks/$IPYNB_FILE"

# echo -e "$txtgrn \n*** STEP 2: $IPYNB_FILE $txtwht"
# papermill \
#     -p data_path  "/data2/kfeng/deeper_guides_dasatinib/dasatinib_primary_01HD/DASA2/cellranger/outs"  \
#     -p h5_name    "$SECONDARY_A_OUT_PATH/seurat/$H5_FILE_NAME" \
#     -p prj_name   "CATATAC_""$PRJ_NAME" \
#     -p protosp_path "$PROTOSPACER_PATH" \
#     "$SECOND_STEP" \
#     "$OUT_PATH"


# # STEP 3 -----------------------------------------------------------------------------
# IPYNB_FILE="03_more_qc_and_viz.R.ipynb"
# THIRD_STEP="$SECONDARY_A_SCRIPT_PATH"/"$IPYNB_FILE"
# OUT_PATH="$SECONDARY_A_OUT_PATH/notebooks/$IPYNB_FILE"

# echo -e "$txtgrn \n*** STEP 3: $IPYNB_FILE $txtwht"
# papermill \
#    -p h5_name    "$SECONDARY_A_OUT_PATH/seurat/$H5_FILE_NAME" \
#    -p protosp_path "$PROTOSPACER_PATH" \
#    "$THIRD_STEP" \
#    "$OUT_PATH"


# # STEP 4 -----------------------------------------------------------------------------
# IPYNB_FILE="04_match_guide2barcode_04.ipynb"
# FOURTH_STEP="$SECONDARY_A_SCRIPT_PATH"/"$IPYNB_FILE"
# OUT_PATH="$SECONDARY_A_OUT_PATH/notebooks/$IPYNB_FILE"

# echo -e "$txtgrn \n*** STEP 4: $IPYNB_FILE $txtwht"
# papermill \
#     -p h5_name    "$SECONDARY_A_OUT_PATH/seurat/$H5_SINGLETS_FILE_NAME" \
#     -p protosp_path "$PROTOSPACER_PATH" \
#     -p library     "$FQ_PREFIX" \
#     -p tsv_path    "$PRIMARY_A_OUT_PATH/tsv" \
#     "$FOURTH_STEP" \
#     "$OUT_PATH"

# STEP 5 -----------------------------------------------------------------------------
 IPYNB_FILE="05_demux_guides_00.ipynb"
 FIFTH_STEP="$SECONDARY_A_SCRIPT_PATH"/"$IPYNB_FILE"
 OUT_PATH="$SECONDARY_A_OUT_PATH/notebooks/$IPYNB_FILE"

 echo -e "$txtgrn \n*** STEP 5: $IPYNB_FILE $txtwht"
 papermill \
     -p h5_name    "$SECONDARY_A_OUT_PATH/seurat/$H5_SINGLETS_FILE_NAME" \
     -p protosp_path "$PROTOSPACER_PATH_05" \
     -p library     "$FQ_PREFIX" \
     -p tsv_path    "$PRIMARY_A_OUT_PATH/tsv" \
     "$FIFTH_STEP" \
     "$OUT_PATH" \
     --kernel ir

