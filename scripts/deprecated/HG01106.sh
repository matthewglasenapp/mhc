#!/bin/bash
#SBATCH --job-name=HG01106
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=HG01106_%j.out
#SBATCH --error=HG01106_%j.err
#SBATCH --mem=60G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1-0
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

# After pbmm2, use samtools to throw away everything not aligned to chromosome 6.  

# 1. pbmarkdup -j <threads> --rmdup <input.bam> <output.fastq>
# 3. cutadapt or fastplong
# 3. pbmm2
# 4. DeepVariant
# 5. pbsv discover
# 6. pbsv call

# Document the size of data (average per sample) at each step.
# One sample takes X amount of time on Y core with Z CPU/RAM 

CCS="/hb/groups/cornejo_lab/HLA_hybrid_capture/Pacbio/20240711_Twist-HLA-Panel/HiFiBam/m84039_240622_113450_s1.hifi_reads.bc1009--bc1009.bam"
REF="/hb/scratch/mglasena/MHC/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
#TRF_BED="human_GRCh38_no_alt_analysis_set.trf.bed"
TRF_BED="test_tandem_repeat.bed"

module load miniconda3
conda activate analysis

# bam2fastq -j 24 $CCS -o HG01106

# pbmarkdup -j 20 --rmdup HG01106.fastq.gz HG01106_dedup.fastq

# gzip HG01106_dedup.fastq

# module load cutadapt/4.4
# cutadapt -j 24 -n 3 -g AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -a "A{10}N{90}" -o HG01106_trimmed_polyA.fastq.gz HG01106_dedup.fastq.gz

# module load miniconda3
# conda activate analysis

# pbmm2 align -j 20 $REF HG01106_trimmed_polyA.fastq.gz HG01106.hg38.bam --sort --log-level INFO --unmapped --bam-index BAI --rg '@RG\tID:m84039_240622_113450_s1\tSM:HG01106'

# samtools view -@ 20 -b HG01106.hg38.bam chr6 > HG01106.chr6.bam
# samtools index HG01106.chr6.bam

pbsv discover --region chr6 --tandem-repeats $TRF_BED HG01106.chr6.bam HG01106.chr6.svsig.gz
tabix -c '#' -s 3 -b 4 -e 4 HG01106.chr6.svsig.gz
pbsv call -j 20 --region chr6 --hifi $REF HG01106.chr6.svsig.gz HG01106.var.vcf

bgzip -c HG01106.var.vcf > HG01106.var.vcf.gz
tabix -p vcf HG01106.var.vcf.gz
bcftools view -r chr6:28000000-34000000 HG01106.var.vcf.gz -o HG01106_MHC.vcf
bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%SVLEN\n' HG01106_MHC.vcf > HG01106_MHC.bed
awk -v OFS="\t" '{ if (NR>1) print $1, $2-1, $3, $4, $5, $6 }' HG01106_MHC.bed > HG01106_MHC_fixed.bed

sort -k1,1 -k2,2n HG01106_MHC_fixed.bed > HG01106_MHC_fixed.sorted.bed

echo -e "CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN\tGene_Chrom\tGene_Start\tGene_End\tGene_ID" > HG01106_SV.tsv
bedtools intersect -a HG01106_MHC_fixed.sorted.bed -b mhc_genes.bed -loj >> HG01106_SV.tsv


