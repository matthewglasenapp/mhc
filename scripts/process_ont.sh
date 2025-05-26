#!/bin/bash
#SBATCH --job-name=process_ont
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=process_ont_%A_%a.out
#SBATCH --error=process_ont_%A_%a.err
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=3-0
#SBATCH --array=0-32
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

echo "=============================="
echo "   ONT Processing Script   "
echo "=============================="

echo -e "\n\n"

echo "Using Python: $(which python3)"
echo "Conda Environment: $(conda info --envs | grep '*')"

echo -e "\n\n"

array_id=$SLURM_ARRAY_TASK_ID
export array_id

module load gatk
module load dorado

python3 -u process_ont.py