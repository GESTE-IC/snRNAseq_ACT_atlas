#!/bin/bash
#
# 10X SpaceRanger pipeline
# version: spaceranger-1.3.1

input_path="/path/to/fastq"
spaceranger="/path/to/space/ranger/spaceranger"
transcriptome_ref="path/to/refdata-gex-GRCh38-2020-A"
probes_ref="/path/to/probe_sets/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv"
ID_sample="sample_to_process"
output_path="path/to/output"

$spaceranger count --id=$ID_sample
                  --transcriptome=$transcriptome_ref \ 
                  --probe-set=$probes_ref \
                  --fastqs=$input_path \
                  --sample=$ID_sample \ 
                  --image="path/to/your/image.tif" \
                  --slide=V10S29-130 \ #slide ID
                  --area=B1 \ #Capture area
                  --loupe-alignment="/my/path/to/loupe/alignment/file/V10S29-130-B1-mysamplename.json" \
                  --localcores=8 \ #allowed cores in localmode
                  --localmem=64 \ #allowed memory in localmode
