#!/bin/bash
set -e
RUN_DIR="/home/vsevim/prj/tf/tf2-novaseq/data/atac_and_guidecap/"

CUR_DIR=$(pwd)
cd $RUN_DIR
mkdir -p fastqc

fastqc ./*.fastq.gz -o fastqc -t 10

cd "$CUR_DIR"
