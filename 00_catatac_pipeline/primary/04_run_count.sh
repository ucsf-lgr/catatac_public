#!/bin/bash
set -euo pipefail
slack() {
    curl -X POST -H 'Content-type: application/json' \
         --data "{\"text\": \"$1\" }" \
         https://hooks.slack.com/services/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
}

CR_ARC_PATH="/home/vsevim/software/cellranger-arc-2.0.1/cellranger-arc"

# cellranger-arc count --id=Condition7 \
#                     --reference=/home/vsevim/prj/refs/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
#                     --libraries=/home/vsevim/prj/tf/tf2-novaseq/workflow/workflow/scripts/library7.csv \
#                     --localcores=64 \
#                     --localmem=750
echo "USAGE: 04_run_count.sh OUTPUT_DIR LIBRARY_TSV [DRY_RUN] " 


REFERENCE=/home/vsevim/prj/refs/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
CONDITION_NAME=$1
LIBRARY_TSV=$2
DRY_RUN=$3

if [[ "$DRY_RUN" == "DRY" ]];  then
    ${CR_ARC_PATH} count --id="$CONDITION_NAME" \
                    --reference="$REFERENCE" \
                    --libraries="$LIBRARY_TSV" \
                    --localcores=64 \
                    --localmem=500 \
                    --dry
else
    ${CR_ARC_PATH} count \
		    --no-bam \
		    --id="$CONDITION_NAME" \
                    --reference="$REFERENCE" \
                    --libraries="$LIBRARY_TSV" \
                    --localcores=64 \
                    --localmem=500 
fi
