#!/bin/bash
#SBATCH --job-name=SV_benchmark
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=SV_benchmark_%j.out
#SBATCH --error=SV_benchmark_%j.err
#SBATCH --mem=25G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1-0
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --account=pi-jkoc

truth_vcf="/hb/scratch/mglasena/MHC/concordance/GIAB_benchmark/HG002_SV/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz"
query_vcf="HG002.dedup.trimmed.hg37.chr6.SV.norm.vcf.gz"
confident_regions="/hb/scratch/mglasena/MHC/concordance/GIAB_benchmark/HG002_SV/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.bed"
regions_file="/hb/scratch/mglasena/test_pacbio/merged_hla_3.bed"
reference_fasta="/hb/scratch/mglasena/test_pacbio/reference/GRCh37/hs37d5.fa"
output_prefix="hap_py_results/HG002"
rtg_path="/hb/home/mglasena/.conda/envs/happy/bin/rtg"
rtg_template="/hb/scratch/mglasena/MHC/concordance/hap_py_input/rtg_sdf_grch37"
threads=8

export HGREF=$reference_fasta

# Don't use -f $confident_regions

hap.py $truth_vcf $query_vcf -R $regions_file -r $reference_fasta -o $output_prefix --engine vcfeval --engine-vcfeval-path $rtg_path --engine-vcfeval-template $rtg_template --threads $threads

#SV Analyzer
# https://github.com/nhansen/SVanalyzer/blob/master/docs/svbenchmark.rst

# truth_vcf="/hb/scratch/mglasena/MHC/concordance/GIAB_benchmark/HG002_SV/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz"
# query_vcf="/hb/scratch/mglasena/test_pacbio/processed_data/GRCh37/pbsv_vcf/HG002.dedup.trimmed.hg37.chr6.SV.vcf.gz"
# reference_fasta="/hb/scratch/mglasena/test_pacbio/reference/GRCh37/hg19.fa"
# svanalyzer benchmark --ref $reference_fasta --test $query_vcf --truth $truth_vcf