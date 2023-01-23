#!/bin/bash

# activate appropriate env
source ~/.bashrc
conda activate CellBender

cd /dfs3b/swaruplab/shared_lab/cross-disorder/count_matrices/AD_Mathys_2019

# get a list of all samples:
sample_dir="/dfs3b/swaruplab/shared_lab/cross-disorder/count_matrices/AD_Mathys_2019/"
samples=("D17-8800" "D17-8754")

for i in "${samples[@]}"
do
    echo $i

    # set input and output files:
    infile=$sample_dir$i/counts_unfiltered/
    outfile=$sample_dir$i/counts_unfiltered/cellbender.h5

    # run cellbender
    cellbender remove-background \
       --input $infile \
       --output $outfile \
       --expected-cells 300 \
       --total-droplets-included 15000 \
       --epochs 150 \
       --fpr 0.01 \
       --cuda

done
