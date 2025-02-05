#!/bin/bash
#SBATCH --job-name=cutadapt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=cutadapt_%j.out
#SBATCH --error=cutadapt_%j.err
#SBATCH --mem=60G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1-0
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

module load miniconda3
conda activate analysis

#BAM="/hb/groups/cornejo_lab/HLA_hybrid_capture/Pacbio/20240711_Twist-HLA-Panel/HiFiBam/m84039_240622_113450_s1.hifi_reads.bc1001--bc1001.bam"
#bam2fastq -j 24 $BAM -o HG002

#module load fastqc/0.12.1
#fastqc --memory 5000 --adapters adapters.txt HG002.fastq.gz

# module load cutadapt/4.4
# cutadapt \
# -j 24 \
# -n 3 \
# --poly-a \
# -g AGATGTGTATAAGAGACAG \
# -a CTGTCTCTTATACACATCT \
# -o HG002_trimmed_polyA.fastq.gz \
# HG002.fastq.gz
# #-a "A{10}N{90}"

# fastplong \
# --thread 24 \
# -Q \
# -L \
# -s AGATGTGTATAAGAGACAG \
# -e CTGTCTCTTATACACATCT \
# --distance_threshold .1 \
# --trim_poly_x \
# --in HG002.fastq.gz \
# --out HG002_trimmed.fastq.gz

module load fastqc/0.12.1
fastqc --memory 5000 HG002_trimmed.fastq.gz