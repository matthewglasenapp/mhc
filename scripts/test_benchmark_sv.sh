# Define reference fasta variable
reference_fasta="/hb/scratch/mglasena/test_pacbio/reference/GRCh37/hs37d5.fa"
tandem_repeat_bed="/hb/scratch/mglasena/test_pacbio/repeats_bed/human_hs37d5.trf.bed"
truth_vcf="/hb/scratch/mglasena/MHC/concordance/GIAB_benchmark/HG002_SV/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz"
regions_file="/hb/scratch/mglasena/test_pacbio/merged_hla.bed"
fastq_file="/hb/scratch/mglasena/test_pacbio/processed_data/fastq_rmdup_cutadapt/HG002.dedup.trimmed.fastq.gz"

# Download truvari
# git clone https://github.com/spiralgenetics/truvari
# (cd truvari; git reset --hard 600b4ed7)


# Map non-ACGT characters to N
# sed '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' "$reference_fasta" > ref/human_hs37d5.fasta

# Download tandem repeat bed annotations
curl -s https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_hs37d5.trf.bed > ref/human_hs37d5.trf.bed

# Map fastq files to new reference
pbmm2 align -j 16 $reference_fasta $fastq_file HG002.dedup.trimmed.hg38.bam --sort --log-level INFO --unmapped --bam-index BAI --rg @RG\tID:m84039_240622_113450_s1\tSM:HG002 

# Filter reads for chromosome 6
samtools view -@ 16 -b HG002.dedup.trimmed.hg38.bam 6 > HG002.dedup.trimmed.hg37.chr6.bam
samtools index HG002.dedup.trimmed.hg37.chr6.bam

# Run pbsv
pbsv discover --region 6 --tandem-repeats $tandem_repeat_bed HG002.dedup.trimmed.hg37.chr6.bam HG002.dedup.trimmed.hg37.chr6.svsig.gz
tabix -c '#' -s 3 -b 4 -e 4 HG002.dedup.trimmed.hg37.chr6.svsig.gz

pbsv call -j 16 --region 6 --hifi $reference_fasta HG002.dedup.trimmed.hg37.chr6.svsig.gz hg2.pbsv.vcf -t INS,DEL
bgzip hg2.pbsv.vcf
tabix hg2.pbsv.vcf.gz

# truvari/truvari.py -f ref/human_hs37d5.fasta -b $truth_vcf --includebed $regions_file --passonly --giabreport -r 1000 -p 0.00 -c hg2.pbsv.vcf.gz -o bench

# --passonly
truvari bench -f ref/human_hs37d5.fasta -b $truth_vcf --includebed mhc_3.bed -r 25000 --chunksize 25000 --sizemin 10 -S 10 -p 0.00 --extend 1000 -c hg2.pbsv.vcf.gz -o bench
