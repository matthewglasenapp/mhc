#!/bin/bash
#SBATCH --job-name=process_pacbio
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=process_pacbio_%A_%a.out
#SBATCH --error=process_pacbio_%A_%a.err
#SBATCH --mem=25G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=3-0
#SBATCH --array=0-32%9
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

array_id=$SLURM_ARRAY_TASK_ID
export array_id

python3 -u process_pacbio.py