#!/bin/bash
#SBATCH --job-name=process_pb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=process_pb_%j.out
#SBATCH --error=process_pb_%j.err
#SBATCH --mem=60G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1-0
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

# module load miniconda3
# conda activate analysis

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

#module load fastqc/0.12.1
#fastqc --memory 5000 HG002_trimmed.fastq.gz

# module load cutadapt/4.4
# cutadapt \
# -j 24 \
# -a "A{10}N{90}" \
# -o HG002_trimmed_polyA.fastq.gz \
# HG002.fastq.gz

# module load fastqc/0.12.1
# fastqc --memory 5000 HG002_trimmed_polyA.fastq.gz

REF="/hb/scratch/mglasena/MHC/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
#FASTQ="HG002_trimmed_polyA.fastq.gz"
#OUTFILE="HG002.hg38.bam"

#pbmm2 align -j 20 $REF $FASTQ $OUTFILE --sort --log-level INFO --unmapped --bam-index BAI --rg '@RG\tID:m84039_240622_113450_s1\tSM:HG002'

BAM="HG002.hg38.bam"
OUTFILE="HG002.hg38.svsig.gz"
OUTFILE2="HG002.var.vcf"
TRF_BED="human_GRCh38_no_alt_analysis_set.trf.bed"
# Tandem Repeat Bed https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
#pbsv discover --region chr6 --tandem-repeats $TRF_BED $BAM $OUTFILE
#tabix -c '#' -s 3 -b 4 -e 4 $OUTFILE
#pbsv call -j 20 --region chr6 --hifi $REF $OUTFILE $OUTFILE2

#bgzip -c HG002.var.vcf > HG002.var.vcf.gz
#tabix -p vcf HG002.var.vcf.gz
#bcftools view -r chr6:28000000-34000000 HG002.var.vcf.gz -o HG002_MHC.vcf

#awk -v OFS="\t" '{ if (NR>1) print $1, $2-1, $3, $4, $5, $6 }' HG002_MHC.bed > HG002_MHC_fixed.bed


# bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%SVLEN\n' HG002_MHC.vcf > HG002_MHC.bed
# awk 'BEGIN{OFS="\t"} {if (NR>1) print $1, $2-1, $3, $4, $5, $6}' HG002_MHC.bed > HG002_MHC_fixed.bed
# echo -e "CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN\tGene_Chrom\tGene_Start\tGene_End\tGene_ID" > MHC_variants_intersect.tsv
# bedtools intersect -a HG002_MHC_fixed.bed -b hla_captured_genes.bed -loj >> MHC_variants_intersect.tsv

# REF="/hb/scratch/mglasena/MHC/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
# BAM="HG002.hg38.bam"
# trgt genotype \
#   --genome $REF \
#   --reads $BAM \
#   --repeats polymorphic_repeats.hg38.bed \
#   --output-prefix HG002 \
#   --karyotype XY \
#   --threads 24 \
#   --preset targeted 

# bcftools concat --allow-overlaps \
# /hb/scratch/mglasena/MHC/genotypes/revio/HG002.hg38.revio.vcf.gz \
# HG002.var.vcf.gz \
# HG002.vcf.gz | \
# grep -v "chrX\|chrY" | \
# grep -v "SVTYPE=BND\|SVTYPE=INV\|SVTYPE=DUP" | \
# bcftools norm -D --fasta-ref $REF | \
# bcftools sort | \
# bgzip > merged.vcf.gz
# tabix merged.vcf.gz

REF="/hb/scratch/mglasena/MHC/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
hiphase \
--threads 16 \
--ignore-read-groups \
--global-realignment-cputime 300 \
--reference $REF \
--bam HG002.hg38.bam \
--output-bam HG002.hg38.haplotag.bam \
--vcf merged.vcf.gz \
--output-vcf HG002.hg38.phased.vcf.gz \
--stats-file HG002.phased.stats \
--blocks-file HG002.phased.blocks \
--summary-file HG002.phased.summary
