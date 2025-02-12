#!/bin/bash
#SBATCH --job-name=process_pacbio
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=process_pacbio_%j.out
#SBATCH --error=process_pacbio_%j.err
#SBATCH --mem=25G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=1-0
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

echo "=============================="
echo "   PacBio Processing Script   "
echo "=============================="

echo -e "\n\n"

echo "Using Singularity: $(which singularity)"
echo "Using Python: $(which python3)"
echo "Conda Environment: $(conda info --envs | grep '*')"

echo -e "\n\n"

python3 -u process_pacbio.py