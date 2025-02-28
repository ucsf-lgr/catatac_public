#!/bin/bash
#RUN_PATH="/mnt/LGR_Nextseq/220215_NB552361_0068_AHMGVWBGXK/"
#OUT_PATH="/home/vsevim/prj/catatac/data/bcl2fq_out/"
#RUN_PATH="/mnt/LGR_Nextseq/220414_NB552361_0074_AHJMVNBGXL"
#OUT_PATH="/home/vsevim/prj/catatac/data/Test2/bcl2fastq/"
#RUN_PATH="/mnt/LGR_Nextseq/220420_NB552361_0077_AHK2T5BGXL/"
#OUT_PATH="/home/vsevim/prj/catatac/data/Test2/bcl2fastq2/"
RUN_PATH="/mnt/LGR_Nextseq/220524_NB552361_0085_AHCNJKAFX3/"
OUT_PATH="/home/vsevim/prj/catatac/data/Test3/cond6/atac_and_guidecap/bcl2fastq"
SAMPLESHEET="/home/vsevim/prj/catatac/scripts/samplesheet_test3_cond6_atac_and_guidecap.csv"
#RUN_PATH="/home/vsevim/prj/catatac/data/Test3/cond6/guide_capture/220526_M06776_0049_000000000-KC6GW/"
#OUT_PATH="/home/vsevim/prj/catatac/data/Test3/cond6/guide_capture/bcl2fastq/"
#SAMPLESHEET="/home/vsevim/prj/catatac/scripts/samplesheet_test3_cond6_guidecap.csv"
# RUN_PATH="/home/vsevim/prj/catatac/data/Test3/cond6/atac/220524_NB552361_0085_AHCNJKAFX3/"
# OUT_PATH="/home/vsevim/prj/catatac/data/Test3/cond6/atac/bcl2fastq/"
# SAMPLESHEET="/home/vsevim/prj/catatac/scripts/samplesheet_test3_cond6_atac.csv"

echo "Using --barcode-mismatches 0"
echo "See if bcl2fastq fails without it." 
echo "I included this line bc two Test2 barcodes were too similar."
echo "For condition6 miseq run, use --use-bases-mask Y70,I8,Y24,Y70 instead of Y70,I8,Y16,Y70"
echo "Also, for contition6 use --barcode-mismatches 2, not 0"

/usr/local/bin/bcl2fastq \
  --use-bases-mask Y70,I8,Y16,Y70 \
  --create-fastq-for-index-reads \
  --minimum-trimmed-read-length=8 \
  --mask-short-adapter-reads=8 \
  --ignore-missing-positions \
  --ignore-missing-controls \
  --ignore-missing-filter \
  --ignore-missing-bcls \
  --barcode-mismatches 0 \
  -r 6 -w 6 \
  -R ${RUN_PATH} \
  --output-dir=${OUT_PATH} \
  --interop-dir=${RUN_PATH}"/InterOp/" \
  --sample-sheet=${SAMPLESHEET} 
  
  #2>"$OUT_PATH/bcl2fq.out"

