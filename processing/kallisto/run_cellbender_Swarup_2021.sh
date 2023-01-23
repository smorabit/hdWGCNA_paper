#!/bin/bash

# activate appropriate env
source ~/.bashrc
conda activate CellBender

cd /dfs3b/swaruplab/shared_lab/cross-disorder/count_matrices/AD_Swarup/

END=18
for index in $(seq 1 $END);
do
    echo $index

    # get the current file to download from the filepaths txt file:
    settings=$(head -n $index /dfs3b/swaruplab/shared_lab/cross-disorder/bin/Swarup_2021/cellbender_settings.csv | tail -n 1)
    cur_sample=$(echo $settings | cut -d ',' -f 1)
    cells=$(echo $settings | cut -d ',' -f 3)
    total_drops=$(echo $settings | cut -d ',' -f 4)

    echo $settings
    echo $cur_sample
    echo $cells
    echo $total_drops


    # set input and output files:
    infile=$cur_sample/counts_unfiltered/
    outfile=$cur_sample/counts_unfiltered/cellbender.h5

    # run cellbender
    cellbender remove-background \
       --input $infile \
       --output $outfile \
       --expected-cells $cells \
       --total-droplets-included $total_drops \
       --epochs 150 \
       --fpr 0.01 \
       --cuda

done

echo $settings
