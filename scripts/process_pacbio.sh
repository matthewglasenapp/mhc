#!/bin/bash
#SBATCH --job-name=process_pacbio
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=process_pacbio_%j.out
#SBATCH --error=process_pacbio_%j.err
#SBATCH --mem=60G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1-0
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

module purge
module load singularity-ce/singularity-ce.4.1.4

module load miniconda3
conda activate analysis

echo "Using Singularity: $(which singularity)"
echo "Using Python: $(which python3)"
echo "Conda Environment: $(conda info --envs | grep '*')"

python3 -u process_pacbio.py