#!/bin/bash
set -euo pipefail

txtred=$'\e[0;31m' # Red
txtgrn=$'\e[0;32m' # Green
txtylw=$'\e[0;33m' # Yellow
txtwht=$'\e[0;37m'

#CONDITION="condition_08"
FQ_PREFIX="CRISPR8"
CONDITION="$FQ_PREFIX"
PRJ_NAME="CATATAC_4_""$CONDITION"
H5_FILE_NAME="$PRJ_NAME"".h5seurat"
H5_SINGLETS_FILE_NAME="$PRJ_NAME""_ONLY_SINGLETS.h5seurat"
SECONDARY_A_SCRIPT_PATH="/home/vsevim/prj/workflows/catatac/secondary/"
PRIMARY_A_OUT_PATH="/home/vsevim/prj/tf/CATATAC_4/analysis/primary/$CONDITION"
SECONDARY_A_OUT_PATH="/home/vsevim/prj/tf/CATATAC_4/analysis/secondary/$CONDITION"
MAIN_DATA_PATH="/home/vsevim/prj/tf/CATATAC_4/data/atac_and_guidecap"
PROTOSPACER_PATH="$PRIMARY_A_OUT_PATH/resources/protospacers_4set.tsv"
PRELIM_RUN=true # used to run part way through 05 to select guides
# after running this script once and editing your ps file, updated the below path and change PRELIM_RUN to false
# this gets added after running 05 once
PROTOSPACER_PATH_05="$PRIMARY_A_OUT_PATH/resources/ps_dasatinib_screen3.tsv"
GUIDE_CALL="YES"

mkdir -p "$SECONDARY_A_OUT_PATH/seurat/"
mkdir -p "$SECONDARY_A_OUT_PATH/notebooks/"



# STEP 1 ----------------------------------------------------------------------------
IPYNB_FILE=01_assign_guides_to_cells.py.ipynb
FIRST_STEP="$SECONDARY_A_SCRIPT_PATH/$IPYNB_FILE"
OUT_PATH="$SECONDARY_A_OUT_PATH/notebooks/$IPYNB_FILE"

echo -e "$txtgrn \n*** STEP 1: $IPYNB_FILE $txtwht"
papermill \
    -p library     "$FQ_PREFIX" \
    -p fastq_path  "$MAIN_DATA_PATH" \
    -p tsv_path    "$PRIMARY_A_OUT_PATH/tsv" \
    "$FIRST_STEP" \
    "$OUT_PATH"


# STEP 2 -----------------------------------------------------------------------------
IPYNB_FILE="02_seurat.R.ipynb"
SECOND_STEP="$SECONDARY_A_SCRIPT_PATH"/"$IPYNB_FILE"
OUT_PATH="$SECONDARY_A_OUT_PATH/notebooks/$IPYNB_FILE"

echo -e "$txtgrn \n*** STEP 2: $IPYNB_FILE $txtwht"
papermill \
    -p data_path  "$PRIMARY_A_OUT_PATH/cellranger/outs"  \
    -p h5_name    "$SECONDARY_A_OUT_PATH/seurat/$H5_FILE_NAME" \
    -p prj_name   "CATATAC_""$PRJ_NAME" \
    -p protosp_path "$PROTOSPACER_PATH" \
    "$SECOND_STEP" \
    "$OUT_PATH"


# STEP 3 -----------------------------------------------------------------------------
IPYNB_FILE="03_more_qc_and_viz.R.ipynb"
THIRD_STEP="$SECONDARY_A_SCRIPT_PATH"/"$IPYNB_FILE"
OUT_PATH="$SECONDARY_A_OUT_PATH/notebooks/$IPYNB_FILE"

echo -e "$txtgrn \n*** STEP 3: $IPYNB_FILE $txtwht"
papermill \
    -p h5_name    "$SECONDARY_A_OUT_PATH/seurat/$H5_FILE_NAME" \
    -p protosp_path "$PROTOSPACER_PATH" \
    "$THIRD_STEP" \
    "$OUT_PATH"


# STEP 4 -----------------------------------------------------------------------------
IPYNB_FILE="04_match_guide2barcode_04.ipynb"
FOURTH_STEP="$SECONDARY_A_SCRIPT_PATH"/"$IPYNB_FILE"
OUT_PATH="$SECONDARY_A_OUT_PATH/notebooks/$IPYNB_FILE"

echo -e "$txtgrn \n*** STEP 4: $IPYNB_FILE $txtwht"
papermill \
    -p h5_name    "$SECONDARY_A_OUT_PATH/seurat/$H5_SINGLETS_FILE_NAME" \
    -p protosp_path "$PROTOSPACER_PATH" \
    -p library     "$FQ_PREFIX" \
    -p tsv_path    "$PRIMARY_A_OUT_PATH/tsv" \
    "$FOURTH_STEP" \
    "$OUT_PATH"

# STEP 5 -----------------------------------------------------------------------------
IPYNB_FILE="05_demux_guides_NGC.ipynb"
FIFTH_STEP="$SECONDARY_A_SCRIPT_PATH"/"$IPYNB_FILE"
OUT_PATH="$SECONDARY_A_OUT_PATH/notebooks/${IPYNB_FILE}"

echo -e "$txtgrn \n*** STEP 5: $IPYNB_FILE $txtwht"
papermill \
    -p h5_name    "$SECONDARY_A_OUT_PATH/seurat/$H5_SINGLETS_FILE_NAME" \
    -p protosp_path "$PROTOSPACER_PATH" \
    -p library     "$FQ_PREFIX" \
    -p tsv_path    "$PRIMARY_A_OUT_PATH/tsv" \
    -p run_guide_caller "$GUIDE_CALL" \
    "$FIFTH_STEP" \
    "$OUT_PATH"
    

