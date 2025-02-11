#!/bin/bash
#SBATCH --job-name=process_pacbio
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=process_pacbio_%A_%a.out
#SBATCH --error=process_pacbio_%A_%a.err
#SBATCH --mem=60G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1-0
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

python3 -u process_pacbio.py