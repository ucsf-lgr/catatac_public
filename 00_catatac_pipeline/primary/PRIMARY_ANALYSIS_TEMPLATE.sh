#!/bin/bash
set -euo pipefail

txtred=$'\e[0;31m' # Red
txtgrn=$'\e[0;32m' # Green
txtylw=$'\e[0;33m' # Yellow
txtwht=$'\e[0;37m'

FQ_PREFIX="CRISPR8"
CONDITION="$FQ_PREFIX"
DRY_RUN="NOT" #"DRY_RUN"
PRIMARY_A_SCRIPT_PATH="/home/vsevim/prj/workflows/catatac/primary"
PRIMARY_A_OUT_PATH="/home/vsevim/prj/tf/CATATAC_4/analysis/primary/$CONDITION/"
FASTQ_PATH="/home/vsevim/prj/tf/CATATAC_4/data/atac_and_guidecap"
SEQUENCER="NOVASEQ"
UMI_LEN=12
ORIG_LIB_CSV_PATH="/home/vsevim/prj/tf/CATATAC_4/resources"
LIB_CSV_NAME="library_cond7_catatac_4.csv"
PS_FILE_NAME="protospacers_4set.tsv"
HAMMING_DISTANCE_CS=0 # KEEP SET TO 0 FOR NOW. mismatch allowance for capture sequences 
HAMMING_DISTANCE=1 # 0 1 OR 2
N_THREADS=64 # num of threads to use for sort

# Create output dirs
mkdir -p $PRIMARY_A_OUT_PATH/cellranger
mkdir -p $PRIMARY_A_OUT_PATH/tsv
mkdir -p $PRIMARY_A_OUT_PATH/resources
cp "$ORIG_LIB_CSV_PATH/$LIB_CSV_NAME" "$PRIMARY_A_OUT_PATH/resources/"
cp "$ORIG_LIB_CSV_PATH/$PS_FILE_NAME" "$PRIMARY_A_OUT_PATH/resources/"

# Run pipeline steps
echo "$txtgrn *** STEP 1 *** $txtwht"
"$PRIMARY_A_SCRIPT_PATH/02_reads2tsv.sh" \
        $FASTQ_PATH \
        $PRIMARY_A_OUT_PATH/tsv \
        $FQ_PREFIX \
        $SEQUENCER \
        $HAMMING_DISTANCE_CS

echo -e "$txtgrn \n*** STEP 2 *** $txtwht"
"$PRIMARY_A_SCRIPT_PATH/03_select_reads_by_guide_with_mismatch.sh" \
        $FQ_PREFIX \
        $PRIMARY_A_OUT_PATH/tsv \
        $UMI_LEN \
        $PRIMARY_A_OUT_PATH/resources \
        $PRIMARY_A_OUT_PATH/resources/$PS_FILE_NAME\
        $HAMMING_DISTANCE \
        $FASTQ_PATH \
        $N_THREADS

echo -e "$txtgrn \n*** STEP 3 *** $txtwht"
cd "$PRIMARY_A_OUT_PATH"
"$PRIMARY_A_SCRIPT_PATH/04_run_count.sh" \
        cellranger \
        "$PRIMARY_A_OUT_PATH/resources/$LIB_CSV_NAME" \
        "$DRY_RUN"

echo -e "$txtgrn \n*** Completed *** $txtwht"
echo -e "$txtylw \n Moving cellranger output to $PRIMARY_A_OUT_PATH/cellranger $txtwht"
mv ./cellranger "$PRIMARY_A_OUT_PATH/cellranger"
        

# mkdir -p $PRIMARY_A_OUT_PATH/cellranger
# mkdir -p $PRIMARY_A_OUT_PATH/tsv
# mkdir -p $PRIMARY_A_OUT_PATH/resources
# cp "$ORIG_LIB_CSV_PATH/$LIB_CSV_NAME" "$PRIMARY_A_OUT_PATH/resources/"
# cp "$ORIG_LIB_CSV_PATH/$PS_FILE_NAME" "$PRIMARY_A_OUT_PATH/resources/"


# echo "$txtgrn *** STEP 1 *** $txtwht"
# "${PRIMARY_A_SCRIPT_PATH}/02_reads2tsv.sh" \
#         $FASTQ_PATH \
#         $PRIMARY_A_OUT_PATH/tsv \
#         $FQ_PREFIX \
#         $SEQUENCER

# echo -e "$txtgrn \n*** STEP 2 *** $txtwht"
# "${PRIMARY_A_SCRIPT_PATH}/03_select_reads_by_guide.sh" \
#         $FQ_PREFIX \
#         $PRIMARY_A_OUT_PATH/tsv \
#         $UMI_LEN \
#         $PRIMARY_A_OUT_PATH/resources \
#         $PRIMARY_A_OUT_PATH/resources/$PS_FILE_NAME

# echo -e "$txtgrn \n*** STEP 3 *** $txtwht"
# cd "$PRIMARY_A_OUT_PATH"
# "${PRIMARY_A_SCRIPT_PATH}/04_run_count.sh" \
#         cellranger \
#         "{$PRIMARY_A_OUT_PATH}/resources/${LIB_CSV_NAME}" \
#         "$DRY_RUN"

# mv ./cellranger "$PRIMARY_A_OUT_PATH/cellranger"
        
