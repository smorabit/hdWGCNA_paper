#!/bin/bash

# activate appropriate env
source ~/.bashrc
conda activate CellBender

cd /dfs3b/swaruplab/shared_lab/cross-disorder/count_matrices/AD_Zhou_2020/

# get a list of all samples:
sample_dir="/dfs3b/swaruplab/shared_lab/cross-disorder/count_matrices/AD_Zhou_2020/"
samples=($(ls $sample_dir))

#
# END=28
# for index in $(seq 1 $END);
# do
#     echo $index
#
#     # get the current file to download from the filepaths txt file:
#     settings=$(head -n $index /dfs3b/swaruplab/shared_lab/cross-disorder/bin/Zhou_2020/cellbender_settings.csv | tail -n 1)
#     cur_sample=$(echo $settings | cut -d ',' -f 1)
#     cells=$(echo $settings | cut -d ',' -f 3)
#     total_drops=$(echo $settings | cut -d ',' -f 4)
#
#     echo $settings
#     echo $cur_sample
#     echo $cells
#     echo $total_drops
#
#
#     # set input and output files:
#     infile=$cur_sample/counts_unfiltered/
#     outfile=$cur_sample/counts_unfiltered/cellbender.h5
#
#     # run cellbender
#     cellbender remove-background \
#        --input $infile \
#        --output $outfile \
#        --expected-cells $cells \
#        --total-droplets-included $total_drops \
#        --epochs 150 \
#        --fpr 0.01 \
#        --cuda
#
# done
#
# echo $settings
#
#
#


# samples that I need to re-do
samples=("TWCC-AD7-AD7")




expected_cells=(4000 5000 3000)

for cur in ${samples[@]};
do
  for cells in ${expected_cells[@]};
  do
      echo "$i $cells"

      mkdir $sample_dir$cur/counts_unfiltered/$cells

      # set input and output files:
      infile=$sample_dir$cur/counts_unfiltered/
      outfile=$sample_dir$cur/counts_unfiltered/$cells/cellbender.h5

      # run cellbender
      cellbender remove-background \
         --input $infile \
         --output $outfile \
         --expected-cells $cells \
         --total-droplets-included 15000 \
         --epochs 150 \
         --fpr 0.01 \
         --cuda

  done
done
