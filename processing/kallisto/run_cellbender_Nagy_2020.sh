#!/bin/bash

# activate appropriate env
source ~/.bashrc
conda activate CellBender

cd /dfs3b/swaruplab/shared_lab/cross-disorder/count_matrices/MDD_Nagy_2020/

# get a list of all samples:
sample_dir="/dfs3b/swaruplab/shared_lab/cross-disorder/count_matrices/MDD_Nagy_2020/"
samples=($(ls $sample_dir))

# these runs are part of the same samples:
# SRR12391930 SRR12391931
# SRR12391932 SRR12391933

END=34
for index in $(seq 1 $END);
do
    echo $index

    # get the current file to download from the filepaths txt file:
    settings=$(head -n $index /dfs3b/swaruplab/shared_lab/cross-disorder/bin/Nagy_2020/cellbender_settings.csv | tail -n 1)
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

# echo $settings


#bad_samples=("SRR12391909" "SRR12391919" "SRR12391932" "SRR12391933")
