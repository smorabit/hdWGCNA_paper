#!/bin/bash
#SBATCH --job-name=velm
#SBATCH -p standard
#SBATCH -A vswarup_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --error=slurm-%J.err
#SBATCH --mem 16G
#SBATCH --array=0-53
#SBATCH --time=72:00:00

source ~/.bashrc
conda activate scrublet

sample_dir="/dfs3b/swaruplab/shared_lab/cross-disorder/count_matrices/ASD_Velmeshev_2019/"
samples=($(ls $sample_dir))

let index="$SLURM_ARRAY_TASK_ID"
sample=${samples[$index]}
echo $sample


indir=$sample_dir$sample"/counts_unfiltered/"
infile=$indir"cellbender_filtered.h5"

# run if we have run cellbender
mkdir $sample_dir$sample"/scanpy_preprocessing_cellbender"
mkdir $sample_dir$sample"/scanpy_preprocessing_cellbender/figures"
mkdir $sample_dir$sample"/scanpy_preprocessing_cellbender/data"

python /dfs3b/swaruplab/shared_lab/cross-disorder/bin/scanpy_preprocessing.py \
  --infile $infile \
  --indir $indir \
  --outdir $sample_dir$sample"/scanpy_preprocessing_cellbender/" \
  --mt 5 \
  --maxcount 20000 \
  --mincount 500 \
  --leidenres 0.5 \
  --npcs 30 \
  --neighbors 20 \
  --mindist 0.25 \
  --velo "True"
