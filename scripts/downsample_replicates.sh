#!/bin/bash
#SBATCH --job-name=downsample
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=downsample_%A_%a.out
#SBATCH --error=downsample_%A_%a.err
#SBATCH --mem=25G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=7-0
#SBATCH --array=0-19
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

array_id=$SLURM_ARRAY_TASK_ID
export array_id

python3 -u downsample_replicates.py