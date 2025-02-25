import subprocess
import shutil
import sys

#==============================
# Download necessary files for analysis

# Download hg19 reference with decoys
# mkdir ref
# curl -s ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz > ref/human_hs37d5.fasta.gz
# gunzip ref/human_hs37d5.fasta.gz

# Map non-ACGT characters to N
# sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' ref/human_hs37d5.fasta

# Download tandem repeat bed annotations
# curl -s https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_hs37d5.trf.bed > ref/human_hs37d5.trf.bed
#==============================

max_threads = 12

reference_fasta = "ref/human_hs37d5.fasta"
tandem_repeat_bed = "ref/human_hs37d5.trf.bed"
truth_vcf = "/hb/scratch/mglasena/MHC/concordance/GIAB_benchmark/HG002_SV/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz"
#regions_file = "mhc_3_hg19.bed"
regions_file = "merged_hla_hg19.bed"
#regions_file = "test.bed"
fastq_file = "/hb/scratch/mglasena/test_pacbio/processed_data/fastq_rmdup_cutadapt/HG002.dedup.trimmed.fastq.gz"
HG002_RG_string = r'"@RG\tID:m84039_240622_113450_s1\tSM:HG002"'

# hg38 HLA Class III region
# MCCD1 start (BED): 31528961
# BTNL2 stop (BED): 32407181
# 6	31528961	32407181

# Lift over hg 19 coordinates from hg38 using CrossMap 
# CrossMap bed hg38ToHg19.over.chain.gz mhc_3_hg38.bed mhc_3_hs37d5.bed
# 6	31496738	32374958

# Or take hg19 coordinates for MCCD1 - BTNL2 directly from ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz
# 6	31496493	32374905

# Ensure all required tools are installed and executable
def check_required_commands():    
	print("Checking the installation status of the required bioinformatics tools!")

	required_commands = [
		"bgzip",
		"pbmm2",
		"pbsv",
		"samtools",
		"tabix",
		"truvari"
	]

	missing_commands = []
	for command in required_commands:
		if shutil.which(command) is None:
			missing_commands.append(command)
	if len(missing_commands) != 0:
		print("Error: Missing the following commands: {}".format(", ".join(missing_commands)))
		sys.exit(1)
	else:
		print("All tools required are installed!")
		print("\n\n")

# Map fastq files to hg19 reference genome
def align_to_reference():
	#pbmm2 align -j 16 $reference_fasta $fastq_file HG002.dedup.trimmed.hg38.bam --sort --log-level INFO --unmapped --bam-index BAI --rg @RG\tID:m84039_240622_113450_s1\tSM:HG002 

	output_bam = "HG002.dedup.trimmed.hg19.bam"

	pbmm2_cmd = "pbmm2 align -j {max_threads} {reference_fasta} --sort --log-level INFO --unmapped --bam-index BAI  {input_file} {output_file} --rg {read_group_string}".format(max_threads = max_threads, reference_fasta = reference_fasta, input_file = fastq_file, output_file = output_bam, read_group_string = HG002_RG_string)

	subprocess.run(pbmm2_cmd, shell=True, check=True)

	# Filter reads for chromosome 6
	#samtools view -@ 16 -b HG002.dedup.trimmed.hg38.bam 6 > HG002.dedup.trimmed.hg37.chr6.bam

	samtools_cmd = "samtools view -@ {max_threads} -b {input_bam} 6 > {output_bam}".format(max_threads = max_threads, input_bam = output_bam, output_bam = "HG002.dedup.trimmed.hg19.chr6.bam")

	subprocess.run(samtools_cmd, shell=True, check=True)

	# samtools index HG002.dedup.trimmed.hg37.chr6.bam
	index_cmd = "samtools index HG002.dedup.trimmed.hg19.chr6.bam"

	subprocess.run(index_cmd, shell=True, check=True)

def run_pbsv():
	# Run pbsv
	# pbsv discover --region 6 --tandem-repeats $tandem_repeat_bed HG002.dedup.trimmed.hg37.chr6.bam HG002.dedup.trimmed.hg37.chr6.svsig.gz
	
	input_bam = "HG002.dedup.trimmed.hg19.chr6.bam"
	output_svsig = "HG002.dedup.trimmed.hg19.chr6.svsig.gz"
	output_vcf = "HG002.dedup.trimmed.hg19.chr6.vcf"

	pbsv_discover_cmd = "pbsv discover --region 6 --tandem-repeats {tandem_repeat_bed} {input_bam} {output_svsig}".format(tandem_repeat_bed = tandem_repeat_bed, input_bam = input_bam, output_svsig = output_svsig)
	
	subprocess.run(pbsv_discover_cmd, shell=True, check=True)

	# tabix -c '#' -s 3 -b 4 -e 4 HG002.dedup.trimmed.hg37.chr6.svsig.gz
	index_cmd = "tabix -c '#' -s 3 -b 4 -e 4 {output_svsig}".format(output_svsig = output_svsig)

	subprocess.run(index_cmd, shell=True, check=True)
	
	# pbsv call -j 16 --region 6 --hifi $reference_fasta HG002.dedup.trimmed.hg37.chr6.svsig.gz hg2.pbsv.vcf -t INS,DEL
	# bgzip hg2.pbsv.vcf
	# tabix hg2.pbsv.vcf.gz
	pbsv_call_cmd = "pbsv call -j {max_threads} --region 6 --hifi {reference_fasta} {input_svsig} {output_vcf}".format(max_threads = max_threads, reference_fasta = reference_fasta, input_svsig = output_svsig, output_vcf = output_vcf)

	subprocess.run(pbsv_call_cmd, shell=True, check=True)

	bgzip_cmd = "bgzip {input_file}".format(input_file = output_vcf)
	tabix_cmd = "tabix {input_file}".format(input_file = output_vcf + ".gz")
	subprocess.run(bgzip_cmd, shell=True, check=True)
	subprocess.run(tabix_cmd, shell=True, check=True)

def run_truvari():
	# --passonly
	#truvari bench -f ref/human_hs37d5.fasta -b $truth_vcf --includebed mhc_3.bed -r 25000 --chunksize 25000 --sizemin 10 -S 10 -p 0.00 --extend 1000 -c hg2.pbsv.vcf.gz -o bench
	
	input_vcf = "HG002.dedup.trimmed.hg19.chr6.vcf.gz"

	# --extend 1000
	truvari_cmd = "truvari bench -f {reference_fasta} -b {truth_vcf} --includebed {bed_file} -r 1000 -p 0 --sizemin 10 --sizefilt 10 -c {input_vcf} -o {output_dir}".format(reference_fasta = reference_fasta, truth_vcf = truth_vcf, bed_file = regions_file, input_vcf = input_vcf, output_dir = "bench_truvari")

	subprocess.run(truvari_cmd, shell=True, check=True)

def run_svanalyzer():
	input_vcf = "HG002.dedup.trimmed.hg19.chr6.vcf.gz"

	# svanalyzer benchmark --ref $reference_fasta --test $query_vcf --truth $truth_vcf
	
	#svanalyzer_cmd = "svanalyzer benchmark --ref {reference_fasta} --test {input_vcf} --truth {truth_vcf}".format(reference_fasta = reference_fasta, input_vcf = input_vcf, truth_vcf = truth_vcf)
	#subprocess.run(svanalyzer_cmd, shell=True, check=True)

def main():
	#check_required_commands()
	#align_to_reference()
	#run_pbsv()
	run_truvari()
	#run_svanalyzer()

if __name__ == "__main__":
	main()

