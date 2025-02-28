#!/bin/bash
txtred=$'\e[0;31m' # Red
txtgrn=$'\e[0;32m' # Green
txtylw=$'\e[0;33m' # Yellow
txtwht=$'\e[0;37m'

# Specify the prefix for the guide capture (CRISPR) fastq
FQ_PREFIX="TF2CRISPR2" 
CONDITION="$FQ_PREFIX"
PRIMARY_A_OUT_PATH="/home/vsevim/prj/tf/CATATAC_4/analysis/primary/$CONDITION/"
FASTQ_PATH="/home/vsevim/prj/tf/CATATAC_4/data/atac_and_guidecap_renamed/"
SEQUENCER="NOVASEQ"
UMI_LEN=12
ORIG_LIB_CSV_PATH="/home/vsevim/prj/tf/CATATAC_4/resources"
LIB_CSV_NAME="library_TF2_x2.csv"
PS_FILE_NAME="protospacers_TF2.csv"
DRY_RUN="NOT" #"DRY_RUN"

mkdir -p $PRIMARY_A_OUT_PATH/cellranger
mkdir -p $PRIMARY_A_OUT_PATH/tsv
mkdir -p $PRIMARY_A_OUT_PATH/resources
cp "$ORIG_LIB_CSV_PATH/$LIB_CSV_NAME" "$PRIMARY_A_OUT_PATH/resources/"
cp "$ORIG_LIB_CSV_PATH/$PS_FILE_NAME" "$PRIMARY_A_OUT_PATH/resources/"

echo "$txtgrn *** STEP 1 *** $txtwht"
~/prj/workflows/catatac/primary/02_reads2tsv.sh \
        "$FASTQ_PATH" \
        "$PRIMARY_A_OUT_PATH/tsv" \
        $FQ_PREFIX \
        $SEQUENCER

echo -e "$txtgrn \n*** STEP 2 *** $txtwht"
~/prj/workflows/catatac/primary/03_select_reads_by_guide.sh \
        "$FQ_PREFIX" \
        "$PRIMARY_A_OUT_PATH/tsv" \
        $UMI_LEN \
        $PRIMARY_A_OUT_PATH/resources \
        $PRIMARY_A_OUT_PATH/resources/$PS_FILE_NAME

~/prj/workflows/catatac/primary/04_run_count.sh \
        cellranger \
        $PRIMARY_A_OUT_PATH/resources/$LIB_CSV_NAME \
        $DRY

mv ./cellranger "$PRIMARY_A_OUT_PATH/cellranger"