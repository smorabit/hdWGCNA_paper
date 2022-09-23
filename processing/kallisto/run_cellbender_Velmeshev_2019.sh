#!/bin/bash

# activate appropriate env
source ~/.bashrc
conda activate CellBender

cd /dfs3b/swaruplab/shared_lab/cross-disorder/count_matrices/ASD_Velmeshev_2019/

# get a list of all samples:
sample_dir="/dfs3b/swaruplab/shared_lab/cross-disorder/count_matrices/ASD_Velmeshev_2019/"
samples=($(ls $sample_dir))


END=55
for index in $(seq 1 $END);
do
    echo $index

    # get the current file to download from the filepaths txt file:
    settings=$(head -n $index /dfs3b/swaruplab/shared_lab/cross-disorder/bin/Velmeshev_2019/cellbender_settings.csv | tail -n 1)
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

# need to re-do with fewer expected cells?
# idk these knee-plots don't look great lmao
# should look at these again to see which I might want to re-do

# stopped looking at sample 389
bad_samples=("SRR9262918" "SRR9262924" "SRR9262926" "SRR9262929" "SRR9262939" "SRR9262946" "SRR9262947" "SRR9262952" "SRR9262953" "SRR9262955" "SRR264388" "SRR264389" "SRR264392")
redo_samples=("SRR9262918" "SRR9262924" "SRR9262926" "SRR9262946" "SRR9262947" "SRR9262955" "SRR264389" )
