#!/bin/bash
#SBATCH --job-name=genotype
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=genotype_%A_%a.out
#SBATCH --error=genotype_%A_%a.err
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=3-0
#SBATCH --array=0-31%5
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

# For running deepvariant
module load bcftools
module load singularity-ce/singularity-ce.4.1.4

# For running Clair3
#module load miniconda3/3.12
#conda activate clair3

array_id=$SLURM_ARRAY_TASK_ID
export array_id

python3 -u genotype.py