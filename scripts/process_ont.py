import os
import subprocess
import shutil
import sys
import time

max_threads = 6

# set input directory to the current working directory where the script should be run
input_dir = os.getcwd()

# Set the output directory to processed_data/
output_dir = os.path.join(input_dir, "processed_data")
os.makedirs(output_dir, exist_ok=True)

# Use reference fasta with no alternate contigs.
reference_fasta = os.path.join(input_dir, "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")

tandem_repeat_bed = os.path.join(input_dir, "human_GRCh38_no_alt_analysis_set.trf.bed")
chr6_bed = os.path.join(input_dir, "chr6.bed")

# Set mapped chr6 reads threshold at which variant calling should not proceed
min_reads_sample = 100

barcode_config = "/hb/scratch/mglasena/test_ont/sample_barcode_arrangement.txt"
barcode_file = "/hb/scratch/mglasena/test_ont/sample_barcodes.fa"
basecalled_reads = "/hb/groups/cornejo_lab/HLA_hybrid_capture/06_25_24_R1041_LIG_Cornejo_EXP26/Cornejo/06_25_24_R1041_LIG_Cornejo_EXP26_1_drd0.7.2_sup5.0.0.bam"

demux_prefix = "b2f9e1ada541ad6c3b470699dfbdd70ff26e092f_"

barcode_sample_dict = {
'MY_CUSTOM_KIT_barcode01': 'HG002', 'MY_CUSTOM_KIT_barcode02': 'HG003', 'MY_CUSTOM_KIT_barcode03': 'HG004', 'MY_CUSTOM_KIT_barcode04': 'HG005', 'MY_CUSTOM_KIT_barcode05': 'NA24694', 'MY_CUSTOM_KIT_barcode06': 'NA24695', 'MY_CUSTOM_KIT_barcode07': 'HG01106', 'MY_CUSTOM_KIT_barcode08': 'HG01258', 'MY_CUSTOM_KIT_barcode09': 'HG01891', 'MY_CUSTOM_KIT_barcode10': 'HG01928', 'MY_CUSTOM_KIT_barcode11': 'HG02055', 'MY_CUSTOM_KIT_barcode12': 'HG02630', 'MY_CUSTOM_KIT_barcode13': 'HG03579', 'MY_CUSTOM_KIT_barcode14': 'NA19240', 'MY_CUSTOM_KIT_barcode15': 'NA20129', 'MY_CUSTOM_KIT_barcode16': 'NA21309', 'MY_CUSTOM_KIT_barcode17': 'HG03492', 'MY_CUSTOM_KIT_barcode18': 'IHW09071', 'MY_CUSTOM_KIT_barcode19': 'IHW09021', 'MY_CUSTOM_KIT_barcode20': 'IHW09175', 'MY_CUSTOM_KIT_barcode21': 'IHW09049', 'MY_CUSTOM_KIT_barcode22': 'IHW09117', 'MY_CUSTOM_KIT_barcode23': 'IHW09118', 'MY_CUSTOM_KIT_barcode24': 'IHW09122', 'MY_CUSTOM_KIT_barcode25': 'IHW09125', 'MY_CUSTOM_KIT_barcode26': 'IHW09251', 'MY_CUSTOM_KIT_barcode27': 'IHW09359', 'MY_CUSTOM_KIT_barcode28': 'IHW09364', 'MY_CUSTOM_KIT_barcode29': 'IHW09245', 'MY_CUSTOM_KIT_barcode30': 'IHW09409', 'MY_CUSTOM_KIT_barcode31': 'IHW09198', 'MY_CUSTOM_KIT_barcode32': 'IHW09200', 'MY_CUSTOM_KIT_barcode33': 'IHW09224'
}

clair3_model_path = "/hb/home/mglasena/.conda/envs/clair3/bin/models/r941_prom_sup_g5014"

# Dictionary of samples to process.
# Sample ID: Read Group String
sample_dict = {
	"HG002" : "@RG\tID:m84039_240622_113450_s1\tSM:HG002",
	"HG003" : "@RG\tID:m84039_240622_113450_s1\tSM:HG003",
	"HG004" : "@RG\tID:m84039_240622_113450_s1\tSM:HG004",
	"HG005" : "@RG\tID:m84039_240622_113450_s1\tSM:HG005",
	"HG01106" : "@RG\tID:m84039_240622_113450_s1\tSM:HG01106",
	"HG01258" : "@RG\tID:m84039_240622_113450_s1\tSM:HG01258",
	"HG01891" : "@RG\tID:m84039_240622_113450_s1\tSM:HG01891",
	"HG01928" : "@RG\tID:m84039_240622_113450_s1\tSM:HG01928",
	"HG02055" : "@RG\tID:m84039_240622_113450_s1\tSM:HG02055",
	"HG02630" : "@RG\tID:m84039_240622_113450_s1\tSM:HG02630",
	"HG03492" : "@RG\tID:m84039_240622_113450_s1\tSM:HG03492",
	"HG03579" : "@RG\tID:m84039_240622_113450_s1\tSM:HG03579",
	"IHW09021" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09021",
	"IHW09049" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09049",
	"IHW09071" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09071",
	"IHW09117" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09117",
	"IHW09118" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09118",
	"IHW09122" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09122",
	"IHW09125" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09125",
	"IHW09175" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09175",
	"IHW09198" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09198",
	"IHW09200" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09200",
	"IHW09224" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09224",
	"IHW09245" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09245",
	"IHW09251" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09251",
	"IHW09359" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09359",
	"IHW09364" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09364",
	"IHW09409" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09409",
	"NA19240" : "@RG\tID:m84039_240622_113450_s1\tSM:NA19240",
	"NA20129" : "@RG\tID:m84039_240622_113450_s1\tSM:NA20129",
	"NA21309" : "@RG\tID:m84039_240622_113450_s1\tSM:NA21309",
	"NA24694" : "@RG\tID:m84039_240622_113450_s1\tSM:NA24694",
	"NA24695" : "@RG\tID:m84039_240622_113450_s1\tSM:NA24695"
}

longphase = "/hb/home/mglasena/software/longphase/longphase_linux-x64"

# Ensure all required tools are installed and executable
def check_required_commands():    
	print("Checking the installation status of the required bioinformatics tools!")

	required_commands = [
		"bcftools",
		"bgzip",
		"dorado",
		"fastqc",
		"gatk",
		"minimap2",
		"pigz",
		"porechop_abi",
		"samtools",
		"sniffles",
		"tabix",
		"whatshap"
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

def run_dorado():
	dorado_cmd = "dorado demux -o {output_dir} --emit-summary --barcode-arrangement {barcode_config} --barcode-sequences {barcode_sequences} --kit-name MY_CUSTOM_KIT --threads {threads} {basecalled_reads}".format(output_dir = output_dir, barcode_config = barcode_config, barcode_sequences = barcode_file, threads = max_threads, basecalled_reads = basecalled_reads)
   
	subprocess.run(dorado_cmd, shell=True, check=True)

def rename_demux_bams():
	raw_bam_dir = os.path.join(output_dir, "raw_bam")
	os.makedirs(raw_bam_dir, exist_ok=True)
	for key, value in barcode_sample_dict.items():
		input_file = os.path.join(input_dir, demux_prefix + key + ".bam")
		output_file = os.path.join(output_dir, "raw_bam", value + ".bam")
		shutil.copy(input_file, output_file)

class Samples:
	raw_bam_dir = os.path.join(output_dir, "raw_bam")
	raw_fastq_dir = os.path.join(output_dir, "raw_fastq")
	fastq_porechop_dir = os.path.join(output_dir, "fastq_porechop")
	mapped_bam_dir = os.path.join(output_dir, "mapped_bam")
	clair3_dir = os.path.join(output_dir, "clair3_vcf")
	sniffles_dir = os.path.join(output_dir, "sniffles")
	whatshap_phased_vcf_dir = os.path.join(output_dir, "phased_vcf_whatshap")
	longphase_phased_vcf_dir = os.path.join(output_dir, "phased_vcf_longphase")

	def __init__(self, sample_ID, read_group_string):
		self.sample_ID = sample_ID
		self.unmapped_bam = os.path.join(Samples.raw_bam_dir, self.sample_ID + ".bam")
		self.read_group_string = read_group_string

		for directory in [Samples.raw_bam_dir, Samples.raw_fastq_dir, Samples.fastq_porechop_dir, Samples.mapped_bam_dir, Samples.clair3_dir, Samples.sniffles_dir, Samples.whatshap_phased_vcf_dir, Samples.longphase_phased_vcf_dir]:
			os.makedirs(directory, exist_ok=True)
		
		print(f"Processing Sample {sample_ID}!")
		print("\n\n")
	
	# Convert BAM file of unmapped HiFi (ccs) reads to FASTQ format for marking duplicates and trimming adapters
	def convert_bam_to_fastq(self):
		print("Converting raw reads to fastq format using gatk SamToFastq!")
		print("SamToFastq input file: {}".format(self.unmapped_bam))

		output_fastq = os.path.join(Samples.raw_fastq_dir, self.sample_ID + ".fastq.gz")

		gatk_SamToFastq_cmd = "gatk SamToFastq -I {input_file} -F {output_file}".format(input_file = self.unmapped_bam, output_file = output_fastq)
		
		subprocess.run(gatk_SamToFastq_cmd, shell=True, check=True)

	def run_porechop_abi(self):
		print("Removing adapters and barcodes with porechop_abi!")
		
		input_fastq = os.path.join(Samples.raw_fastq_dir, self.sample_ID + ".fastq.gz")
		
		print("porechop_abi input file: {}".format(input_fastq))

		porechop_threads = 4

		output_fastq = os.path.join(Samples.fastq_porechop_dir, self.sample_ID + ".porechop.fastq")
		
		porechop_cmd = "porechop_abi --ab_initio -i {input_file} -t {threads} -o {output_file} --format fastq".format(input_file = input_fastq, threads = porechop_threads, output_file = output_fastq)

		subprocess.run(porechop_cmd, shell=True, check=True)

	def trim_reads(self):
		print("Trimiming reads with ProwlerTrimmer!")

		prowler_trimmer = "/hb/home/mglasena/software/ProwlerTrimmer/TrimmerLarge.py"
		input_fastq_dir = Samples.fastq_porechop_dir + "/"
		input_fastq_file = self.sample_ID + ".porechop.fastq"

		prowler_trimmer_cmd = 'python3 {prowler_trimmer} -i {input_dir} -f {input_file} -o {output_dir} -m "D" -q 20'.format(prowler_trimmer = prowler_trimmer, input_dir = input_fastq_dir, input_file = input_fastq_file, output_dir = input_fastq_dir)

		subprocess.run(prowler_trimmer_cmd, shell=True, check=True)

		trimmed_fastq = os.path.join(Samples.fastq_porechop_dir, self.sample_ID + ".porechopTrimLT-U0-D20W100L100R0.fastq")
		pigz_cmd = "pigz -f -p {threads} {input_file}".format(threads = max_threads, input_file = trimmed_fastq)
		subprocess.run(pigz_cmd, shell=True, check=True)

	def align_to_reference(self):
		print("Aligning reads to GRCh38 reference genome with minimap2!")
		
		input_fastq = os.path.join(Samples.fastq_porechop_dir, self.sample_ID + ".porechopTrimLT-U0-D20W100L100R0.fastq.gz")

		print("minimap2 input file: {}".format(input_fastq))
		
		output_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.bam")

		# Use minimap2 instead of pbmm2
		minimap_threads = int(max_threads * 2 / 3)
		samtools_threads = max_threads - minimap_threads
		minimap_rg_string = "'{}'".format(self.read_group_string.replace("\t", "\\t"))

		minimap2_cmd = "minimap2 -t {minimap_threads} -ax map-ont {reference_genome} {input_file} -R {rg_string} | samtools sort -@ {samtools_threads} -o {output_file}".format(minimap_threads = minimap_threads, reference_genome = reference_fasta, input_file = input_fastq, rg_string = minimap_rg_string, samtools_threads = samtools_threads, output_file = output_bam)

		index_bam = "samtools index {input_file}".format(input_file = output_bam)
		
		subprocess.run(minimap2_cmd, shell=True, check=True)
		subprocess.run(index_bam, shell=True, check=True)

		print("Mapped bam written to: {}".format(output_bam))
		print("\n\n")

	def mark_duplicates(self):
		input_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.bam")
		output_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.mrkdup.bam")
		metrics_file = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.mrkdup.metrics.txt")
		temp_dir = os.path.join(output_dir, "mark_duplicates")
		os.makedirs(temp_dir, exist_ok=True)
		
		mark_duplicates_cmd = "gatk MarkDuplicates -I {input_file} -O {output_file} --TMP_DIR {temp_dir} -M {metrics_file} --CREATE_INDEX true".format(input_file = input_bam, output_file = output_bam, temp_dir = temp_dir, metrics_file = metrics_file)
		
		subprocess.run(mark_duplicates_cmd, shell=True, check=True)

	def filter_reads(self):
		print("Excluding BAM records that don't map to chromosome 6!")
		
		input_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.mrkdup.bam")

		print("Samtools input file: {}".format(input_bam))

		output_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.bam")

		# Exclude duplicates and secondary and supplementary alignments
		samtools_cmd = "samtools view -@ {threads} -F 3328 -b {input_file} chr6 > '{output_file}'".format(threads = max_threads, input_file = input_bam, output_file = output_bam)
		
		subprocess.run(samtools_cmd, shell=True, check=True)

		index_cmd = "samtools index {input_file}".format(input_file = output_bam)

		subprocess.run(index_cmd, shell=True, check=True)

		count_reads_cmd = "samtools view -c {input_file}".format(input_file = output_bam)

		read_count = int(subprocess.check_output(count_reads_cmd, shell=True).strip())
		
		print("Filtered BAM records written to: {}".format(output_bam))
		print("\n\n")

		print(read_count)
		return read_count

	def call_variants(self):
		print("Calling SNVs and small INDELS with Clair3!")

		input_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.bam")
		output_dir = os.path.join(Samples.clair3_dir, self.sample_ID)
		os.makedirs(output_dir, exist_ok=True)

		clair3_cmd = "run_clair3.sh --bam_fn={input_file} --ref_fn={reference_genome} --platform=ont --model_path={clair3_model} --output={output_dir} --threads={threads} --sample_name={sample} --bed_fn={bed_file}".format(input_file = input_bam, reference_genome = reference_fasta, clair3_model = clair3_model_path, output_dir = output_dir, threads = max_threads, sample = self.sample_ID, bed_file = chr6_bed)
		
		subprocess.run(clair3_cmd, shell=True, check=True)

		print("VCF written to {}".format(os.path.join(Samples.clair3_dir, self.sample_ID, "merge_output.vcf.gz")))
		print("\n\n")

	def call_structural_variants(self):
		print("Calling structural variants with Sniffles!")
		
		input_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.bam")
		output_vcf = os.path.join(Samples.sniffles_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.sniffles.vcf.gz")

		sniffles_cmd = "sniffles --allow-overwrite -t {threads} --regions {bed_file} -i {input_bam} -v {output_vcf} --tandem-repeats {tandem_repeat_bed}".format(threads = max_threads, bed_file = chr6_bed, input_bam = input_bam, output_vcf = output_vcf, tandem_repeat_bed = tandem_repeat_bed)
		subprocess.run(sniffles_cmd, shell=True, check=True)

	def phase_genotypes_whatshap(self):
		print("Phasing Genotypes with WhatsHap!")

		input_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.bam")
		input_vcf = os.path.join(Samples.clair3_dir, self.sample_ID, "merge_output.vcf.gz")

		print("Input BAM: {}".format(input_bam))
		print("Input VCF: {}".format(input_vcf))

		haplotagged_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.whatshap.haplotag.bam")
		phased_vcf = os.path.join(Samples.whatshap_phased_vcf_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.phased.vcf.gz")
		output_blocks_file = os.path.join(Samples.whatshap_phased_vcf_dir, self.sample_ID + ".phased.haploblocks.txt")
		output_gtf_file = os.path.join(Samples.whatshap_phased_vcf_dir, self.sample_ID + ".phased.haploblocks.gtf")

		whatshap_phase_cmd = "whatshap phase --ignore-read-groups --output {output_file} --reference {reference_genome} {input_vcf} {input_bam}".format(output_file = phased_vcf, reference_genome = reference_fasta, input_vcf = input_vcf, input_bam = input_bam)

		index_cmd = "bcftools index {input_file}".format(input_file = phased_vcf)
		tabix_cmd = "tabix {input_file}".format(input_file = phased_vcf)

		whatshap_haplotag_cmd = "whatshap haplotag --ignore-read-groups --output {output_file} --reference {reference_genome} {input_vcf} {input_bam}".format(output_file = haplotagged_bam, reference_genome = reference_fasta, input_vcf = phased_vcf, input_bam = input_bam)

		whatshap_stats_cmd = "whatshap stats --block-list={block_list_file} --gtf={gtf_file} {input_vcf}".format(block_list_file = output_blocks_file, gtf_file = output_gtf_file, input_vcf = phased_vcf)

		# Log WhatsHap in own output file so it doesn't clog up STDOUT
		whatshap_log = os.path.join(Samples.whatshap_phased_vcf_dir, self.sample_ID + ".whatshap.log")

		with open(whatshap_log, "w") as log_file:
			log_file.write("\n==== Running WhatsHap Phase ====\n")
			subprocess.run(whatshap_phase_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
			subprocess.run(index_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
			subprocess.run(tabix_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
			log_file.write("\n==== Running WhatsHap Haplotag ====\n")
			subprocess.run(whatshap_haplotag_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
			log_file.write("\n==== Running WhatsHap Stats ====\n")
			subprocess.run(whatshap_stats_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

		print("WhatsHap phased VCF written to: {}".format(phased_vcf))
		print("WhatsHap haplotagged BAM written to: {}".format(haplotagged_bam))
		print("WhatsHap phase block gtf written to: {}".format(output_gtf_file))
		print("WhatsHap phase blocks written to: {}".format(output_blocks_file))
		print("\n\n")

	def phase_genotypes_longphase(self):
		print("Phasing Genotypes with LongPhase!")

		input_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.bam")
		input_SNV_vcf = os.path.join(Samples.clair3_dir, self.sample_ID, "merge_output.vcf.gz")
		input_SV_vcf = os.path.join(Samples.sniffles_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.sniffles.vcf.gz")

		print("Input BAM: {}".format(input_bam))
		print("Input SNV VCF: {}".format(input_SNV_vcf))
		print("Input SV VCF: {}".format(input_SV_vcf))

		haplotagged_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.longphase.haplotag.bam")
		phased_vcf = os.path.join(Samples.longphase_phased_vcf_dir, self.sample_ID + ".porechop.trimmed.hg38.rmdup.chr6.longphase.vcf.gz")
		output_blocks_file = os.path.join(Samples.longphase_phased_vcf_dir, self.sample_ID + ".phased.haploblocks.txt")
		output_gtf_file = os.path.join(Samples.longphase_phased_vcf_dir, self.sample_ID + ".phased.haploblocks.gtf")

		phased_vcf_prefix = phased_vcf.split(".vcf.gz")[0]
		longphase_phase_cmd = "{longphase} phase -s {input_snv_vcf} --sv {input_sv_vcf} -b {input_bam} -r {reference_genome} -t {threads} -o {phased_prefix} --ont".format(longphase = longphase, input_snv_vcf = input_SNV_vcf, input_sv_vcf = input_SV_vcf, input_bam = input_bam, reference_genome = reference_fasta, threads = max_threads, phased_prefix = phased_vcf_prefix)

		compress_cmd = "bgzip -f {prefix}.vcf".format(prefix=phased_vcf_prefix)
		index_cmd = "bcftools index {input_file}".format(input_file = phased_vcf)
		tabix_cmd = "tabix {input_file}".format(input_file = phased_vcf)

		haplotagged_bam_prefix = haplotagged_bam.split(".bam")[0]
		longphase_haplotag_cmd = "{longphase} haplotag -r {reference_genome} -s {input_snv_vcf} --sv-file {input_sv_vcf} -b {input_bam} -t {threads} -o {prefix}".format(longphase = longphase, reference_genome = reference_fasta, input_snv_vcf = phased_vcf, input_sv_vcf = input_SV_vcf, input_bam = input_bam, threads = max_threads, prefix = haplotagged_bam_prefix)

		whatshap_stats_cmd = "whatshap stats --block-list={block_list_file} --gtf={gtf_file} {input_vcf}".format(block_list_file = output_blocks_file, gtf_file = output_gtf_file, input_vcf = phased_vcf)

		# Log WhatsHap in own output file so it doesn't clog up STDOUT
		longphase_log = os.path.join(Samples.whatshap_phased_vcf_dir, self.sample_ID + ".longphase.log")

		with open(longphase_log, "w") as log_file:
			log_file.write("\n==== Running LongPhase Phase ====\n")
			subprocess.run(longphase_phase_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
			subprocess.run(compress_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
			subprocess.run(index_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
			subprocess.run(tabix_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
			log_file.write("\n==== Running LongPhase Haplotag ====\n")
			subprocess.run(longphase_haplotag_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)
			log_file.write("\n==== Running WhatsHap Stats ====\n")
			subprocess.run(whatshap_stats_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

		print("LongPhase phased VCF written to: {}".format(phased_vcf))
		print("LongPhase haplotagged BAM written to: {}".format(haplotagged_bam))
		print("LongPhase phase block gtf written to: {}".format(output_gtf_file))
		print("LongPhase phase blocks written to: {}".format(output_blocks_file))
		print("\n\n")

def main():
	# Check that all required tools are installed
	# check_required_commands()

	# run_dorado()
	# rename_demux_bams()

	array_id = os.environ["array_id"]
	print("Array ID: {}".format(array_id))
	sample_ID = list(sample_dict)[int(array_id)]
	sample_read_group_string = sample_dict[sample_ID]
	print(sample_ID)
	start_time = time.time()
	sample = Samples(sample_ID, sample_read_group_string)
	# sample.convert_bam_to_fastq()
	# sample.run_porechop_abi()
	# sample.trim_reads()
	# sample.align_to_reference()
	# sample.mark_duplicates()

	chr6_reads = sample.filter_reads()

	if chr6_reads > min_reads_sample:
		# sample.call_variants()
		# sample.call_structural_variants()
		# sample.phase_genotypes_whatshap()
		sample.phase_genotypes_longphase()
		end_time = time.time()
		elapsed_time = end_time - start_time
		minutes, seconds = divmod(elapsed_time,60)
		print("Processed sampled in {}:{:.2f}!".format(int(minutes), seconds))

	else:
		print("Insufficient reads for variant calling")
		print("Sample {sample_id} had {num_reads} reads!".format(sample_id = sample_ID, num_reads = chr6_reads))

if __name__ == "__main__":
	main()
