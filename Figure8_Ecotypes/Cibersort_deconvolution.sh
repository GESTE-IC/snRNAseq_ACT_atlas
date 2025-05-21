#!/bin/bash

input_path="path/to/cibersort/ref_matrices/"
output_path="path/to/cibersort/results/"
user="user.name"
token_cibersort="token_cibersort"

cd $input_path

ref="refmatrix_celltype_deconv.txt"
dir_name=$(echo $ref | rev | cut -c5- | rev)

query="pseudobulk_celltype_deconv.txt"
output_dir=$output_path"/pseudobulk_celltype_deconv"
sudo docker run -v $input_path:/src/data -v $output_dir:/src/outdir cibersortx/fractions --username $user --token $token_cibersort --single_cell TRUE --refsample $ref --mixture $query --fraction 0.1 --rmbatchSmode TRUE --verbose TRUE
