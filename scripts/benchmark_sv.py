import subprocess
import shutil
import sys
import os

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
regions_dir = "regions/"
#regions_file = "regions/merged_hla_hg19.bed"
#regions_file = "regions/mhc_3_hg19.bed"
#regions_file = "regions/mhc_3_hg19_captured.bed"
regions_file = "regions/50X_hs37d5.bed"
fastq_file = "/hb/scratch/mglasena/test_pacbio/processed_data/fastq_rmdup_cutadapt/HG002.dedup.trimmed.fastq.gz"
HG002_RG_string = r'"@RG\tID:m84039_240622_113450_s1\tSM:HG002"'

# hg38 HLA Class III region
# MCCD1 start (BED): 31528961
# BTNL2 stop (BED): 32407181
# 6	31528961	32407181

# Lift over hg 19 coordinates from hg38 using CrossMap 
# CrossMap bed hg38ToHg19.over.chain.gz mhc_3_hg38.bed mhc_3_hs37d5.bed
# awk '{ $1="6"; print }' OFS="\t" input.bed > output.bed
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
		"svanalyzer",
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
	output_bam = "mapped_bam/HG002.dedup.trimmed.hg19.bam"

	pbmm2_cmd = "pbmm2 align -j {max_threads} {reference_fasta} --sort --log-level INFO --unmapped --bam-index BAI  {input_file} {output_file} --rg {read_group_string}".format(max_threads = max_threads, reference_fasta = reference_fasta, input_file = fastq_file, output_file = output_bam, read_group_string = HG002_RG_string)

	subprocess.run(pbmm2_cmd, shell=True, check=True)

	# Filter reads for chromosome 6
	filtered_bam = "mapped_bam/HG002.dedup.trimmed.hg19.chr6.bam"
	samtools_cmd = "samtools view -@ {max_threads} -b {input_bam} 6 > {output_bam}".format(max_threads = max_threads, input_bam = output_bam, output_bam = filtered_bam)

	subprocess.run(samtools_cmd, shell=True, check=True)

	index_cmd = "samtools index {input_bam}".format(input_bam = filtered_bam)

	subprocess.run(index_cmd, shell=True, check=True)

def run_pbsv():
	# Run pbsv	
	input_bam = "mapped_bam/HG002.dedup.trimmed.hg19.chr6.bam"
	output_svsig = "pbsv_vcf/HG002.dedup.trimmed.hg19.chr6.svsig.gz"
	output_vcf = "pbsv_vcf/HG002.dedup.trimmed.hg19.chr6.vcf"

	pbsv_discover_cmd = "pbsv discover --region 6 --tandem-repeats {tandem_repeat_bed} {input_bam} {output_svsig}".format(tandem_repeat_bed = tandem_repeat_bed, input_bam = input_bam, output_svsig = output_svsig)
	
	subprocess.run(pbsv_discover_cmd, shell=True, check=True)

	index_cmd = "tabix -c '#' -s 3 -b 4 -e 4 {output_svsig}".format(output_svsig = output_svsig)

	subprocess.run(index_cmd, shell=True, check=True)
	
	pbsv_call_cmd = "pbsv call -j {max_threads} --region 6 --hifi {reference_fasta} {input_svsig} {output_vcf}".format(max_threads = max_threads, reference_fasta = reference_fasta, input_svsig = output_svsig, output_vcf = output_vcf)

	subprocess.run(pbsv_call_cmd, shell=True, check=True)

	bgzip_cmd = "bgzip {input_file}".format(input_file = output_vcf)
	tabix_cmd = "tabix {input_file}".format(input_file = output_vcf + ".gz")
	subprocess.run(bgzip_cmd, shell=True, check=True)
	subprocess.run(tabix_cmd, shell=True, check=True)

def run_truvari():	
	input_vcf = "pbsv_vcf/HG002.dedup.trimmed.hg19.chr6.vcf.gz"

	truvari_cmd = "truvari bench -f {reference_fasta} -b {truth_vcf} --includebed {bed_file} --extend 1000 -r 1000 -p 0 --sizemin 10 --sizefilt 10 -c {input_vcf} -o {output_dir}".format(reference_fasta = reference_fasta, truth_vcf = truth_vcf, bed_file = regions_file, input_vcf = input_vcf, output_dir = "truvari")

	subprocess.run(truvari_cmd, shell=True, check=True)

def run_svanalyzer():
	input_vcf = os.path.join(os.getcwd(), "pbsv_vcf/HG002.dedup.trimmed.hg19.chr6.vcf.gz")
	ref = os.path.abspath(reference_fasta)
	bed = os.path.abspath(regions_file)
	os.chdir("svanalyzer/")
	
	svanalyzer_cmd = "svanalyzer benchmark --ref {reference_fasta} --test {input_vcf} --truth {truth_vcf} -minsize 10 --includebed {bed_file} -prefix svanalyzer".format(reference_fasta = ref, input_vcf = input_vcf, truth_vcf = truth_vcf, bed_file = bed)
	subprocess.run(svanalyzer_cmd, shell=True, check=True)

def main():
	check_required_commands()
	
	output_dirs = ["mapped_bam/", "pbsv_vcf/", "svanalyzer/"]
	#for dir in output_dirs:
		#os.makedirs(dir, exist_ok=True)

	#align_to_reference()
	#run_pbsv()
	run_truvari()
	run_svanalyzer()

if __name__ == "__main__":
	main()

