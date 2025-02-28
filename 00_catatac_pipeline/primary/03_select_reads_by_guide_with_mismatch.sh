#!/bin/bash
set -euo pipefail
#IFS=$'\n\t'

echo "USAGE: 03_select_reads_by_guide.sh PREFIX TSV_DIR UMI_LEN RESOURCE_DIR PROTOSPACER_FILE HAMMING_DISTANCE FASTQ_PATH N_THREADS"

# Inputs specified here
LIBRARY=$1   # ("CRISPR6")
TSV_DIR=$2   # "/home/vsevim/prj/catatac/analysis/Test3/cond6/guide_capture/barcode2guide/"
UMI_LEN=$3   # 8 or 12
RESOURCE_DIR=$4
PROTOSPACER_FILE=$5
HAMMING_DISTANCE=$6 # 0, 1, or 2 
FASTQ_PATH=$7 # new for mismatch
N_THREADS=$8 # n_threads for sorting
# .....................
echo $LIBRARY
echo $TSV_DIR
echo $UMI_LEN
HEADER="atac_barcode\tumi\tprotospacer\tguide\tgene\tread_name\tbarcode"


atac_wlist_path="/home/vsevim/software/cellranger-arc-2.0.0/lib/python/atac/barcodes/737K-arc-v1.txt.gz"
gex_wlist_path="/home/vsevim/software/cellranger-arc-2.0.0/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz"

########################################
##### Mismatch vars and functions ######
########################################
# make this a user input 
# mismatch variables, might be better to user specify
CS1_REFERENCE="/home/sfederman/CAP2001/230426_A01831_0041_AH7VW3DSX7_KS01/UMI_count/1mm/reference/CS1.230217_A01102_0560_BHW2WLDSX5_ACLGR01_ATAC.fa"
READ_WITH_CS1="R1"

PS_REFERENCE="/data2/kev/angela_catatac/primary/CRISPR9_01HD/resources/ps_ref.fa"
READ_WITH_PS="R3"
UMI_LENGTH=12
##################################################
TMP_FOLDER="tmp"
OUTPUT_DIR="${LIBRARY}_0${HAMMING_DISTANCE}HD"
##### End of normal params to change ######

# mismatch functions
make_flat_file () {
    # Make a flat file consisting of
    # header    ID of PS/CS hit for this read   UMI

    local REPAIRED_FASTQ_R1="$1"
    local REPAIRED_FASTQ_R2="$2"
    local FLAT_FILE_TABLE="$3"
    local UMI_LENGTH="$4"

    local BASE1=$(basename "$REPAIRED_FASTQ_R1" .fastq)
    local BASE2=$(basename "$REPAIRED_FASTQ_R2" .fastq)

    local FILE1_HEADERS="${TMP_FOLDER}/${BASE1}.TEMP.headers"
    local FILE2_HEADERS="${TMP_FOLDER}/${BASE2}.TEMP.headers"
    local UMI_SEQUENCE="${TMP_FOLDER}/${BASE1}.TEMP.UMI"

    # Make UMI list
    grep -A1 --no-group-separator "^@" "$REPAIRED_FASTQ_R1" | grep -v "^@" | cut -c1-"$UMI_LENGTH" > "$UMI_SEQUENCE"

    grep "^@" "$REPAIRED_FASTQ_R1" | awk '{print $1,$3}' > "$FILE1_HEADERS"
    grep "^@" "$REPAIRED_FASTQ_R2" | awk '{print $3}'    > "$FILE2_HEADERS"

    paste "$FILE1_HEADERS" "$FILE2_HEADERS" "$UMI_SEQUENCE" > "$FLAT_FILE_TABLE"

    #rm "$FILE1_HEADERS"
    #rm "$FILE2_HEADERS"
    #rm "$UMI_SEQUENCE"

}

map_read () {
    local FASTQ="$1"
    local OUTMATCHED="$2"
    local REF="$3"
    local KMER="$4"
    local HAMMING_DISTANCE="$5"
    local STATS="$6"
    local OUTMASKED="$7"
    local LOG="$8"
    #local OUTPUT_DIR="$9"

    # consider adding these options:

    # rename=f            Rename reads to indicate which sequences they matched.
    # findbestmatch=f     (fbm) If multiple matches, associate read with sequence
    #      sharing most kmers.  Reduces speed.
    # rcomp=f             Default is to [t], to look for reverse complement, might need to revcomp the guides/PS and use f to prevent incorrect hits

    bbduk.sh in="$FASTQ" outm="$OUTMATCHED" ref="$REF" k="$KMER" hdist="$HAMMING_DISTANCE" stats="$STATS" rename=t rcomp=f maskmiddle=f 1>"$LOG" 2>&1
    bbduk.sh in="$OUTMATCHED" out="$OUTMASKED" ref="$REF" k="$KMER" hdist="$HAMMING_DISTANCE" kmask=lc rename=t rcomp=f maskmiddle=f 1>>"$LOG" 2>&1
}


if [[ ! -d "$OUTPUT_DIR" ]]
then
    mkdir "$OUTPUT_DIR"
fi

if [[ ! -d "$TMP_FOLDER" ]]
then
    mkdir "$TMP_FOLDER"
fi

#########################################
#### Mismatch vars and functions end ####
#########################################

# Produce atac to gex barcode lookup table
if test -f "$RESOURCE_DIR/atac2gex.txt"; then
    echo "atac2gex.txt found."
else
    echo "atac2gex.txt not found. Generating..."
    cp $atac_wlist_path .
    gunzip 737*.gz
    mv 737K-arc-v1.txt $RESOURCE_DIR/atac.txt

    cp $gex_wlist_path .
    gunzip 737*.gz
    mv 737K-arc-v1.txt $RESOURCE_DIR/gex.txt

    paste -d"\t" $RESOURCE_DIR/atac.txt $RESOURCE_DIR/gex.txt > $RESOURCE_DIR/atac2gex.txt
    rm $RESOURCE_DIR/atac.txt $RESOURCE_DIR/gex.txt
    echo "Done."  
fi

# Read protospacer file
while read targetgene guidename protospacer_seq; do 
   PROTOSPACERS_GENENAMES+=($targetgene)
   PROTOSPACERS_GUIDENAMES+=($guidename)
   PROTOSPACERS+=($protospacer_seq)
done < <(grep -Pv "^#" ${PROTOSPACER_FILE})

# Label reads with protospacer and guide names
for TSV in "$TSV_DIR"/$LIBRARY*_CS.TSV.gz; do
    echo "TSV: $TSV"
    TSV_NAME=$(basename "$TSV" .w_CS.TSV.gz)
    OUT_NAME="$TSV_DIR/$TSV_NAME"".structured.TSV"
    echo "OUTNAME: $OUT_NAME"
    # run mismatch functions
    for R1_FASTQ in "${FASTQ_PATH}"/*"${LIBRARY}"*"${READ_WITH_CS1}"*.fastq.gz
    do
        R2_FASTQ="${R1_FASTQ/${READ_WITH_CS1}/${READ_WITH_PS}}"

        echo "R1_FASTQ (CS): ${R1_FASTQ}"
        echo "R2_FASTQ (PS): ${R2_FASTQ}"

        BASE_R1=$(basename "$R1_FASTQ" .fastq.gz)
        MATCHED_FASTQ_R1="${OUTPUT_DIR}/${BASE_R1}_D${HAMMING_DISTANCE}.matched.fastq"
        MASKED_FASTQ_R1="${OUTPUT_DIR}/${BASE_R1}_D${HAMMING_DISTANCE}.masked.fastq"

        BASE_R2=$(basename "$R2_FASTQ" .fastq.gz)
        MATCHED_FASTQ_R2="${OUTPUT_DIR}/${BASE_R2}_D${HAMMING_DISTANCE}.matched.fastq"
        MASKED_FASTQ_R2="${OUTPUT_DIR}/${BASE_R2}_D${HAMMING_DISTANCE}.masked.fastq"


        STATS_R1="${OUTPUT_DIR}/${BASE_R1}_D${HAMMING_DISTANCE}.stats.txt"
        REFSTATS_R1="${OUTPUT_DIR}/${BASE_R1}_D${HAMMING_DISTANCE}.refstats.txt"

        STATS_R2="${OUTPUT_DIR}/${BASE_R2}_D${HAMMING_DISTANCE}.stats.txt"
        REFSTATS_R2="${OUTPUT_DIR}/${BASE_R2}_D${HAMMING_DISTANCE}.refstats.txt"

        LOG_R1="${OUTPUT_DIR}/${BASE_R1}_D${HAMMING_DISTANCE}.bbduk.log"
        LOG_R2="${OUTPUT_DIR}/${BASE_R2}_D${HAMMING_DISTANCE}.bbduk.log"

        KMER_R1=22
        KMER_R2=20

        #Search for CS1 in R1
        map_read "$R1_FASTQ" "$MATCHED_FASTQ_R1" "$CS1_REFERENCE" "$KMER_R1" "$HAMMING_DISTANCE" "$STATS_R1" "$MASKED_FASTQ_R1" "$LOG_R1" &

        #Search for PS in R2
        map_read "$R2_FASTQ" "$MATCHED_FASTQ_R2" "$PS_REFERENCE" "$KMER_R2" "$HAMMING_DISTANCE" "$STATS_R2" "$MASKED_FASTQ_R2" "$LOG_R2" &
    done

    # Wait for above to complete before merging
    for job in $(jobs -p)
    do
        wait $job
    done

    # Merge separate analysis from R1/R2 in order to find read pairs that match both criteria

    for R1_FASTQ in "$FASTQ_PATH"*R1*.fastq.gz
    do
        R2_FASTQ="${R1_FASTQ/R1/${READ_WITH_PS}}"

        BASE_R1=$(basename "$R1_FASTQ" .fastq.gz)
        MATCHED_FASTQ_R1="${OUTPUT_DIR}/${BASE_R1}_D${HAMMING_DISTANCE}.matched.fastq"
        REPAIRED_FASTQ_R1="${OUTPUT_DIR}/${BASE_R1}_D${HAMMING_DISTANCE}.repaired.fastq"

        BASE_R2=$(basename "$R2_FASTQ" .fastq.gz)
        MATCHED_FASTQ_R2="${OUTPUT_DIR}/${BASE_R2}_D${HAMMING_DISTANCE}.matched.fastq"
        REPAIRED_FASTQ_R2="${OUTPUT_DIR}/${BASE_R2}_D${HAMMING_DISTANCE}.repaired.fastq"

        R2_GUIDECOUNT="${OUTPUT_DIR}/${BASE_R2}_D${HAMMING_DISTANCE}.guidecount.txt"


        SINGLETON="${OUTPUT_DIR}/${BASE_R1}.singletons.fastq"

        #UMI_COUNT="${OUTPUT_DIR}/${BASE_R1}_D${HAMMING_DISTANCE}.UMIcount.fastq"
        LOG_REPAIR="${OUTPUT_DIR}/${BASE_R1}_D${HAMMING_DISTANCE}.repair.log"

        FLAT_FILE_TABLE="${OUTPUT_DIR}/${BASE_R1}_D${HAMMING_DISTANCE}.merged.headers"

        if [[ -s "$MATCHED_FASTQ_R1" ]] && [[ -s "$MATCHED_FASTQ_R2" ]]
        then
            repair.sh in="$MATCHED_FASTQ_R1" in2="$MATCHED_FASTQ_R2" out="$REPAIRED_FASTQ_R1" out2="$REPAIRED_FASTQ_R2" outs="$SINGLETON" repair 1>"$LOG_REPAIR" 2>&1
            grep "^@" "$REPAIRED_FASTQ_R2" | awk '{print $3}' | sort | uniq -c | sort -nr | awk '{print $2"\t"$1}' | sort -k1 > "$R2_GUIDECOUNT"
            make_flat_file "$REPAIRED_FASTQ_R1" "$REPAIRED_FASTQ_R2" "$FLAT_FILE_TABLE" "$UMI_LENGTH"
            #grep "^[ACGTN]" "$REPAIRED_FASTQ_R1" | cut -c1-"$UMI_LENGTH" | sort | uniq -c | sort -nr > "$UMI_COUNT"
        fi
    done
done


echo "-------------------------------------------"
echo "Make summary TSV with bbduk outputs"
# make summary file from bbmap outputs 
# define vars for making all tsv file 
TSVGZ=$(find "$TSV_DIR" -type f -not -name '*w_CS*' -print -quit) # output of Volkan 02 script

# repaired R1 and R2 bbduk outputs
R1=( "$(find "$OUTPUT_DIR" -name '*R1*repaired.fastq')" )
R2=( "$(find "$OUTPUT_DIR" -name '*R3*repaired.fastq')" )

echo "TSVGZ: ${TSVGZ}"
echo "Repaired R1: ${R1}"
echo "Repaired R2: ${R2}"

# header for output file: atac_barcode  umi protospacer guide   gene    read_name   barcode
# umi, from bbduk R1 repaired, first 8 or 12 bases
awk -v len="$UMI_LEN" 'NR%4 == 2 {print substr($0, 1, len)}' "$R1" > "${LIBRARY}"_umi.tmp
echo "Made umi.tmp"
# guide, from bbduk R2 repaired
grep "^@" "$R2" | awk -F' ' '{split($3, a, "="); print a[1]}' > "${LIBRARY}"_guide.tmp
echo "Made guide.tmp"
# protospacer using guide and ps reference
awk '{print $2,$3}' "$PROTOSPACER_FILE" > "${LIBRARY}"_ps_list.tmp
awk 'NR==FNR {ps[$1]=$2; next} {print $0"\t"ps[$1]}' "${LIBRARY}"_ps_list.tmp "${LIBRARY}"_guide.tmp | awk '{print $2}' > "${LIBRARY}"_ps.tmp
echo "Made ps.tmp"
# gene, capture string from guide before _
grep "^@" "$R2" | awk -F' ' '{split($3, a, "_"); print a[1]}' > "${LIBRARY}"_gene.tmp
echo "Made gene.tmp"
# read_name, R1/R2 bbduk repaired
grep "^@" "$R2" | awk -F' ' '{print substr($1, 2)}' > "${LIBRARY}"_read_name.tmp
echo "Made read_name.tmp"
# combine
paste "${LIBRARY}"_umi.tmp "${LIBRARY}"_ps.tmp "${LIBRARY}"_guide.tmp "${LIBRARY}"_gene.tmp "${LIBRARY}"_read_name.tmp > "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp1
echo "Made tmp1"
# now sort this by read_name
sort -S 50% --parallel=64 -t $'\t' -k5 "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp1 > "${LIBRARY}"_ALL_STRUCTURED.TSV.sorted.tmp1
echo "tmp1 sorted"
# sort $TSVGZ, slow
zcat "$TSVGZ" | awk -F'\t' '{print $1, $4}' | sort -S 50% --parallel="${N_THREADS}" -t $'\t' -k1 > "${LIBRARY}"_TSVGZ_SORTED.tmp
echo "tsvgz sorted"
# atac_barcode
awk -F'\t' '{print $5}' "${LIBRARY}"_ALL_STRUCTURED.TSV.sorted.tmp1 > "${LIBRARY}"_read_name.sorted.tmp
rm "${LIBRARY}"_read_name.tmp
awk 'NR==FNR {atac[$1]=$2; next} {print $0"\t"atac[$1]}' "${LIBRARY}"_TSVGZ_SORTED.tmp "${LIBRARY}"_read_name.sorted.tmp | awk '{print $2}' > "${LIBRARY}"_atac_barcode.tmp 
echo "Made atac_barcode.tmp"

# add atac_barcode to first column
paste "${LIBRARY}"_atac_barcode.tmp "${LIBRARY}"_ALL_STRUCTURED.TSV.sorted.tmp1 > "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp2
# head "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp2

# add barcode to last column
awk 'NR==FNR {gex[$1]=$2; next} {print $0"\t"gex[$1]}' "${RESOURCE_DIR}/atac2gex.txt" "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp2 > "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp3
# echo ""${LIBRARY}"_ALL_STRUCTURED.TSV.tmp3"
# head "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp3
# wc -l "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp3

# sort by guide
sort -s -k4,4 "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp3 > "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp4
# head "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp4

# add header
echo -e "$HEADER" | cat - "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp4 > "${TSV_DIR}"/"${LIBRARY}"_ALL_STRUCTURED.TSV

# remove temps
rm "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp1 "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp2 "${LIBRARY}"_ALL_STRUCTURED.TSV.sorted.tmp1 \
"${LIBRARY}"_ALL_STRUCTURED.TSV.tmp3 "${LIBRARY}"_ALL_STRUCTURED.TSV.tmp4 "${LIBRARY}"_atac_barcode.tmp "${LIBRARY}"_umi.tmp \
"${LIBRARY}"_ps.tmp "${LIBRARY}"_guide.tmp "${LIBRARY}"_gene.tmp "${LIBRARY}"_read_name.sorted.tmp "${LIBRARY}"_TSVGZ_SORTED.tmp \
"${LIBRARY}"_ps_list.tmp

ls -l -h "$TSV_DIR"/"${LIBRARY}"_ALL_STRUCTURED.TSV
echo Written "$TSV_DIR"/"${LIBRARY}"_ALL_STRUCTURED.TSV

# write counts file
cut -f4 "$TSV_DIR"/"${LIBRARY}"_ALL_STRUCTURED.TSV | \
    fgrep -v guide | sort | uniq -c > \
    "$TSV_DIR"/"${LIBRARY}"_COUNTS.TXT

echo "Counts are in " "$TSV_DIR"/"${LIBRARY}"_COUNTS.TXT
echo "Completed."
