#!/bin/bash
#
# 10X Cellranger pipeline
# version: cellranger-3.1.0

input_path="/path/to/fastq"
cellranger="path/to/cellranger/cellranger-3.1.0/cellranger"
transcriptome_ref="path/to/GRCh38-3.0.0.premrna"
ID_sample="sample_to_process"
output_path="path/to/output"

cd $output_path
$cellranger count --id=$ID_sample \
                  --transcriptome=$transcriptome_ref \
                  --fastqs=$input_path
