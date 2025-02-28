#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

echo "USAGE: ./02_reads2tsv.sh \ "
echo "        FASTQ_DIRECTORY \ "
echo "        OUTPUT_PATH \ "
echo "        PREFIX \ "
echo "        SEQUENCER \ "
echo
printf "This version looks for two capture sequences!!!"
echo Run using \"./02_reads2tsv.sh miseq\" when processing a MiSeq run.

RUN_DIR=$1
OUT_PATH=$2
PREFIX=$3
INSTRUMENT=$4
HD_CS=$5 # make sure this is set to 0 since this hasn't been fully tested yet

echo ${RUN_DIR}
echo ${OUT_PATH}
echo ${INSTRUMENT}
echo

CS1="CAAGTTGATAACGGACTAGCC"
CS2="CAAGTTGTAAACGGACTAGCC"


INSTRUMENT=$(echo $INSTRUMENT | tr 'a-z' 'A-Z')
REMOVE_INTERMEDIATE_FILES=0

revcomp() {
    echo "$1" | tr ACGTacgt TGCAtgca | rev
}

REVCOMP_CS1=$(revcomp "$CS1")

#RUN_NAME=$(basename "$RUN_DIR")
#echo "Run name: $RUN_NAME"

for R1_FASTQ_PATH in "$RUN_DIR"/"$PREFIX"*R1*.fastq.gz; do
    R1_FASTQ_FILENAME=$(basename "$R1_FASTQ_PATH")
    R1_BASE=$(basename "$R1_FASTQ_PATH" .fastq.gz)
    R2_FASTQ_PATH="${R1_FASTQ_PATH//_R1_/_R3_}"
    R2_FASTQ_FILENAME=$(basename "$R2_FASTQ_PATH")
    R2_BASE=$(basename "$R2_FASTQ_PATH" .fastq.gz)
    BARCODE_FASTQ_PATH="${R1_FASTQ_PATH//_R1_/_R2_}"
    ALL_READS_TSV_PATH="$OUT_PATH/${R1_BASE}.TSV"
    ALL_READS_TSV_W_CS_PATH="$OUT_PATH/${R1_BASE}.w_CS.TSV"

    echo "Processing $R1_BASE"
    bioawk -c fastx '{print $name,$seq}' "$R1_FASTQ_PATH" >R1.dummy
    echo "Completed R1"

    bioawk -c fastx '{print $seq}' "$R2_FASTQ_PATH" >R2.dummy
    echo "Completed R2"

    if [[ ${INSTRUMENT} = "MISEQ" ]]; then
        echo "Instrument is $INSTRUMENT"
        echo "Trimming 24 base barcodes..."
        bioawk -c fastx '{print substr($seq,1,16)}' "$BARCODE_FASTQ_PATH" >BARCODE.dummy
    elif [[ ${INSTRUMENT} = "NOVASEQ" ]]; then
        echo "Instrument is $INSTRUMENT"
        echo "Trimming 24 base barcodes..."
        bioawk -c fastx '{print revcomp(substr($seq,9,24))}' "$BARCODE_FASTQ_PATH" >BARCODE.dummy
    else
        bioawk -c fastx '{print revcomp($seq)}' "$BARCODE_FASTQ_PATH" >BARCODE.dummy   
    fi
    echo "Completed R3 (barcode)"

    # 1. Concatenate R1 R2 and Barcode into one file
    # 1.1 Remove read if R1 has no Capture Sequence
    paste R1.dummy R2.dummy BARCODE.dummy >"$ALL_READS_TSV_PATH"
    #awk -v v="$CS1" '$2 ~ v ' "$ALL_READS_TSV_PATH" >"$ALL_READS_TSV_W_CS_PATH"
    if (("$HD_CS" == 0)); then
        # if HD_CS is 0 then we perfect match CS's
        grep "$CS1\|$CS2" "$ALL_READS_TSV_PATH" >"$ALL_READS_TSV_W_CS_PATH"
    else
        # Function to generate patterns with a given HD_CS
        echo "CD hamming distance: ${HD_CS}"
        generate_patterns_with_hamming_distance() {
            local pattern="$1"
            local distance=$2
            local len=${#pattern}
            local patterns=()
            for ((i = 0; i < len; i++)); do
                for c in A C G T; do
                    if [[ "${pattern:i:1}" != "$c" ]]; then
                        new_pattern="${pattern:0:i}$c${pattern:i+1}"
                        patterns+=("$new_pattern")
                        if (( distance > 1 )); then
                            generate_patterns_with_hamming_distance "$new_pattern" $((distance - 1))
                        fi
                    fi
                done
            done
            echo "${patterns[@]}"
        }

        # Generate patterns for CS1 and CS2 with the specified Hamming distance
        CS1_patterns=($(generate_patterns_with_hamming_distance "$CS1" "$HD_CS"))
        CS2_patterns=($(generate_patterns_with_hamming_distance "$CS2" "$HD_CS"))

        # Combine both sets of patterns
        all_patterns=("${CS1_patterns[@]}" "${CS2_patterns[@]}")

        # Use grep with the generated patterns to search for matching lines in the TSV file
        grep -E "$(printf '%s|' "${all_patterns[@]}")" "$ALL_READS_TSV_PATH" > "$ALL_READS_TSV_W_CS_PATH"
    fi
    echo "All complete"
    rm R1.dummy R2.dummy BARCODE.dummy
done

pigz $OUT_PATH/*_L???_R?_*.TSV
