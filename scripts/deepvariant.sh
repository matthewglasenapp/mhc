#!/bin/bash
#SBATCH --job-name=deepvariant
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=deepvariant_%A_%a.out
#SBATCH --error=deepvariant_%A_%a.err
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=3-0
#SBATCH --array=0-31%5
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

module load bcftools
module load singularity-ce/singularity-ce.4.1.4
#module load miniconda3/3.12
#conda activate alignment

array_id=$SLURM_ARRAY_TASK_ID
export array_id

python3 -u deepvariant.py