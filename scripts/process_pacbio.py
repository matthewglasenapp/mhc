import os
import shutil
import subprocess
import sys
import time

"""
Work Flow
1. Check that all required bioinformatics tools are installed and executable
2. Convert raw HiFi reads to FASTQ format bam2fastq
3. Remove PCR duplicates with pbmarkdup
4. Run fastqc on the raw deduplicated reads 
5. Trim adapter sequences and polyA tails with cutadapt
6. Run fastqc on the trimmed reads to verify removal of TE and polyA sequences 
7. Map reads to GRCh38 reference genome with pbmm2
8. Filter for reads mapping to chr6 
9. Call SNVs and small INDELS with DeepVariant
10. Call structural variants with pbsv discover and pbsv call
11. Genotype tandem repeats with pbtrgt
12. Merge vcfs from DeepVariant, pbsv, and pbtrgt
13. Phase the merged vcf with HiPhase
"""

"""
Script Details:

	1. The hg 38 fasta reference "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa" and appropriate index file(s) are contained in /current_working_dir/reference.

	2. A DeepVariant Singularity Image Format (sif) file is contained in /current_working_dir/deepvariant_sif. You must build this file on your own machine. Modify deepvariant_cmd in the call_variants() function if you are not using Singularity to run DeepVariant. 

	3. The raw HiFi reads are contained in /current_working_dir/raw_hifi_reads/. The files are named <sample_ID>.hifi_reads.bam. To customize this, edit the constructor of the Samples class. 

	4. The repeat bed files required for pbsv (human_GRCh38_no_alt_analysis_set.trf.bed) and pbtrgt (polymorphic_repeats.hg38.bed) are contained within /current_working_dir/repeats_bed/

	5. The default is set to use 6 threads. Change the max_threads variable to customize

"""

# set input directory to the current working directory where the script should be run
input_dir = os.getcwd()

# Set the output directory to processed_data/
output_dir = os.path.join(input_dir, "processed_data")
os.makedirs(output_dir, exist_ok=True)

# Input file paths

# Use reference fasta with no alternate contigs.
reference_fasta = os.path.join(input_dir, "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")

# DeepVariant sif file
deepvariant_sif = os.path.join(input_dir, "deepvariant_sif/deepvariant.sif")

# GRCh38 tandem repeat mask file for pbsv
# Downloaded from https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
tandem_repeat_bed = os.path.join(input_dir, "repeats_bed/human_GRCh38_no_alt_analysis_set.trf.bed")

# GRCh38 tandem repeat definition file for pbtrgt
# Downloaded from https://zenodo.org/records/8329210
pbtrgt_repeat_file = os.path.join(input_dir, "repeats_bed/polymorphic_repeats.hg38.bed")

# Dictionary of samples to process.
# Sample ID: Read Group String
sample_dict = {
	"HG002" : "@RG\tID:m84039_240622_113450_s1\tSM:HG002",
	# "HG003" : "@RG\tID:m84039_240622_113450_s1\tSM:HG003",
	# "HG004" : "@RG\tID:m84039_240622_113450_s1\tSM:HG004",
	# "HG005" : "@RG\tID:m84039_240622_113450_s1\tSM:HG005",
	# "HG01106" : "@RG\tID:m84039_240622_113450_s1\tSM:HG01106",
	# "HG01258" : "@RG\tID:m84039_240622_113450_s1\tSM:HG01258",
	# "HG01891" : "@RG\tID:m84039_240622_113450_s1\tSM:HG01891",
	# "HG01928" : "@RG\tID:m84039_240622_113450_s1\tSM:HG01928",
	# "HG02055" : "@RG\tID:m84039_240622_113450_s1\tSM:HG02055",
	# "HG02630" : "@RG\tID:m84039_240622_113450_s1\tSM:HG02630",
	# "HG03492" : "@RG\tID:m84039_240622_113450_s1\tSM:HG03492",
	# "HG03579" : "@RG\tID:m84039_240622_113450_s1\tSM:HG03579",
	# "IHW09021" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09021",
	# "IHW09049" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09049",
	# "IHW09071" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09071",
	# "IHW09117" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09117",
	# "IHW09118" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09118",
	# "IHW09122" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09122",
	# "IHW09125" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09125",
	# "IHW09175" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09175",
	# "IHW09198" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09198",
	# "IHW09200" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09200",
	# "IHW09224" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09224",
	# "IHW09245" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09245",
	# "IHW09251" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09251",
	# "IHW09359" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09359",
	# "IHW09364" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09364",
	# "IHW09409" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09409",
	# "NA19240" : "@RG\tID:m84039_240622_113450_s1\tSM:NA19240",
	# "NA20129" : "@RG\tID:m84039_240622_113450_s1\tSM:NA20129",
	# "NA21309" : "@RG\tID:m84039_240622_113450_s1\tSM:NA21309",
	# "NA24694" : "@RG\tID:m84039_240622_113450_s1\tSM:NA24694",
	# "NA24695" : "@RG\tID:m84039_240622_113450_s1\tSM:NA24695"
}

# Max threads available for parallelization
max_threads = 6

# Minimum reads per sample
# DeepVariant is stalling and not exiting for samples with very few BAM records (e.g., HG01891: 35 mapped reads to chr6)
# Set mapped chr6 reads threshold at which variant calling should not proceed
min_reads_sample = 100000

# Transposase mosaic end binding sequence
# The TE sequence (and its reverse complement) introduced during tagmentation still needs to be removed
# Adapters and barcodes were removed by PacBio with lima
me = "AGATGTGTATAAGAGACAG"
me_rc = "CTGTCTCTTATACACATCT"

# Ensure all required tools are installed and executable
def check_required_commands():    
	print("Checking the installation status of the required bioinformatics tools!")

	required_commands = [
		"bam2fastq",
		"bcftools",
		"bgzip",
		"cutadapt",
		"fastqc",
		"hiphase",
		"pbmarkdup",
		"pbmm2",
		"pbsv",
		"pigz",
		"samtools",
		"singularity",
		"tabix",
		"trgt",
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

class Samples:
	fastq_raw_dir = os.path.join(output_dir, "fastq_raw")
	fastq_rmdup_dir = os.path.join(output_dir, "fastq_rmdup")
	fastq_rmdup_cutadapt_dir = os.path.join(output_dir, "fastq_rmdup_cutadapt")
	mapped_bam_dir = os.path.join(output_dir, "mapped_bam")
	deepvariant_dir = os.path.join(output_dir, "deepvariant_vcf")
	pbsv_dir = os.path.join(output_dir, "pbsv_vcf")
	pbtrgt_dir = os.path.join(output_dir, "pbtrgt_vcf")
	merged_vcf_dir = os.path.join(output_dir, "merged_vcf")
	hiphase_phased_vcf_dir = os.path.join(output_dir, "phased_vcf_hiphase")
	whatshap_phased_vcf_dir = os.path.join(output_dir, "phased_vcf_whatshap")

	def __init__(self, sample_ID, read_group_string):
		self.sample_ID = sample_ID
		self.unmapped_bam = os.path.join(input_dir, "raw_hifi_reads", self.sample_ID + ".hifi_reads.bam")
		self.read_group_string = read_group_string

		for directory in [Samples.fastq_raw_dir, Samples.fastq_rmdup_dir, Samples.fastq_rmdup_cutadapt_dir, Samples.mapped_bam_dir, Samples.deepvariant_dir, Samples.pbsv_dir, Samples.pbtrgt_dir, Samples.merged_vcf_dir, Samples.hiphase_phased_vcf_dir, Samples.whatshap_phased_vcf_dir]:
			os.makedirs(directory, exist_ok=True)
		
		print(f"Processing Sample {sample_ID}!")
		print("\n\n")

	# Convert BAM file of unmapped HiFi (ccs) reads to FASTQ format for marking duplicates and trimming adapters
	def convert_bam_to_fastq(self):
		print("Converting HiFi ccs reads to fastq format using pbtk bam2fastq!")
		print("bam2fastq input file: {}".format(self.unmapped_bam))
		
		os.chdir(Samples.fastq_raw_dir)
		
		bam2fastq_cmd = "bam2fastq -j {threads} {input_file} -o {output_prefix}".format(threads = max_threads, input_file = self.unmapped_bam, output_prefix = self.sample_ID)
		
		subprocess.run(bam2fastq_cmd, shell=True, check=True)
				
		print("Raw fastq reads written to: {}".format(os.path.join(Samples.fastq_raw_dir, self.sample_ID + ".fastq.gz")))
		print("\n\n")

	# Mark PCR duplicates with pbmarkdup
	def mark_duplicates(self):
		print("Removing PCR duplicates using pbmarkdup!")
		
		input_fastq = os.path.join(Samples.fastq_raw_dir, self.sample_ID + ".fastq.gz")
		
		print("pbmarkdup input file: {}".format(input_fastq))

		output_fastq = os.path.join(Samples.fastq_rmdup_dir, self.sample_ID + ".dedup.fastq")

		pbmarkdup_cmd = "pbmarkdup -j {threads} --rmdup {input_file} {output_file}".format(threads = max_threads, input_file = input_fastq, output_file = output_fastq)
		
		subprocess.run(pbmarkdup_cmd, shell=True, check=True)

		gzip_cmd = "pigz -p 8 {input_file}".format(input_file = output_fastq)
		subprocess.run(gzip_cmd, shell=True, check=True)
		
		print("De-duplicated reads written to: {}!".format(output_fastq))
		print("\n\n")

	# Run fastqc on fastq file
	def run_fastqc(self, fastqc_file):
		print("Running fastqc on {}".format(fastqc_file))
		
		fastqc_cmd = "fastqc --memory {memory} {input_fastq}".format(memory = 5000, input_fastq = fastqc_file)
		
		subprocess.run(fastqc_cmd, shell=True, check=True)
		
		print("\n\n")

	# Trim Adapter Sequences and polyA tails
	# Customize this command based on the fastqc output
	def trim_adapters(self):
		print("Trimming adapter sequences with cutadapt!")
		
		input_file = os.path.join(Samples.fastq_rmdup_dir, self.sample_ID + ".dedup.fastq.gz")
		
		print("cutadapt input file: {}".format(input_file))

		output_file = os.path.join(Samples.fastq_rmdup_cutadapt_dir, self.sample_ID + ".dedup.trimmed.fastq.gz")
		
		# The following argument is used to trim polyA tails: -a 'A{{10}}N{{90}}. The double brackets are for python string interpolation
		# me and me_rc are the transposase mosaic end binding sequences defined at the top of the script
		cutadapt_cmd = "cutadapt -j {threads} --quiet -n {num_cuts_allowed} -g {five_prime_adapter} -a {three_prime_adapter} -a 'A{{10}}N{{90}}' -o {output_fastq} {input_fastq}".format(threads = max_threads, num_cuts_allowed = 3, five_prime_adapter = me, three_prime_adapter = me_rc, output_fastq = output_file, input_fastq = input_file)

		subprocess.run(cutadapt_cmd, shell=True, check=True)
		
		print("Deduplicated trimmed reads written to: {}".format(output_file))
		print("\n\n")

	# Align to GRCh38 reference genome with pbmm2
	def align_to_reference(self):
		print("Aligning reads to GRCh38 reference genome with pbmm2!")
		
		input_fastq = os.path.join(Samples.fastq_rmdup_cutadapt_dir, self.sample_ID + ".dedup.trimmed.fastq.gz")

		print("pbmm2 input file: {}".format(input_fastq))
		
		output_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.bam")

		pbmm2_cmd = "pbmm2 align -j {threads} {reference_genome} {input_file} {output_file} --sort --log-level INFO --unmapped --bam-index BAI --rg '{rg_string}'".format(threads = max_threads, reference_genome = reference_fasta, input_file = input_fastq, output_file = output_bam, rg_string = self.read_group_string)

		subprocess.run(pbmm2_cmd, shell=True, check=True)
		
		print("Mapped bam written to: {}".format(output_bam))
		print("\n\n")

	# Filer reads that did not map to chromosome 6
	def filter_reads(self):
		print("Excluding BAM records that don't map to chromosome 6!")
		
		input_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.bam")

		print("Samtools input file: {}".format(input_bam))

		output_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")

		samtools_cmd = "samtools view -@ {threads} -b {input_file} chr6 > '{output_file}'".format(threads = max_threads, input_file = input_bam, output_file = output_bam)
		
		subprocess.run(samtools_cmd, shell=True, check=True)

		index_cmd = "samtools index {input_file}".format(input_file = output_bam)

		subprocess.run(index_cmd, shell=True, check=True)

		count_reads_cmd = "samtools view -c {input_file}".format(input_file = output_bam)

		read_count = int(subprocess.check_output(count_reads_cmd, shell=True).strip())
		
		print("Filtered BAM records written to: {}".format(output_bam))
		print("\n\n")

		return read_count

	# Call SNV with DeepVariant
	def call_variants(self):
		print("Calling SNVs and small INDELS with DeepVariant!")

		input_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")

		print("DeepVariant input file: {}".format(input_bam))
		
		output_vcf = os.path.join(Samples.deepvariant_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.SNV.vcf.gz")
		output_gvcf = os.path.join(Samples.deepvariant_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.SNV.g.vcf.gz")

		bind_paths = [
			f"{Samples.deepvariant_dir}:/data",
			f"{Samples.mapped_bam_dir}:/input",
			f"{os.path.dirname(reference_fasta)}:/reference"
		]

		bind_flags = " ".join("--bind {}".format(path) for path in bind_paths)

		deepvariant_cmd = """
			singularity exec {binds} {sif} /opt/deepvariant/bin/run_deepvariant \
				--model_type=PACBIO \
				--ref=/reference/{ref_filename} \
				--reads=/input/{sample}.dedup.trimmed.hg38.chr6.bam \
				--output_vcf=/data/{sample}.dedup.trimmed.hg38.chr6.SNV.vcf.gz \
				--output_gvcf=/data/{sample}.dedup.trimmed.hg38.chr6.SNV.g.vcf.gz \
				--regions chr6 \
				--num_shards=8
			""".format(
				binds=bind_flags,
				sif=deepvariant_sif,
				ref_filename=os.path.basename(reference_fasta),
				sample=self.sample_ID
				)

		# Log DeepVariant in own output file so it doesn't clog up STDOUT
		deepvariant_log = os.path.join(Samples.deepvariant_dir, self.sample_ID + ".deepvariant.log")

		with open(deepvariant_log, "w") as log_file:
			subprocess.run(deepvariant_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

		print("VCF written to {}".format(output_vcf))
		print("GVCF written to {}".format(output_gvcf))
		print("\n\n")

	# Run pbsv to call structural variants (SV)
	def call_structural_variants(self):
		print("Calling structural variants with pbsv!")

		input_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")

		print("pbsv input file: {}".format(input_bam))
		
		output_svsig = os.path.join(Samples.pbsv_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.svsig.gz")
		output_vcf = os.path.join(Samples.pbsv_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.SV.vcf")

		pbsv_discover_cmd = "pbsv discover --region chr6 --tandem-repeats {tandem_repeat_file} {input_file} {output_file}".format(tandem_repeat_file = tandem_repeat_bed, input_file = input_bam, output_file = output_svsig)
		
		subprocess.run(pbsv_discover_cmd, shell=True, check=True)

		index_svsig_cmd = "tabix -c '#' -s 3 -b 4 -e 4 {svsig_file}".format(svsig_file = output_svsig)
		
		subprocess.run(index_svsig_cmd, shell=True, check=True)

		pbsv_call_cmd = "pbsv call -j {threads} --region chr6 --hifi {reference_genome} {input_file} {output_file}".format(threads = max_threads, reference_genome = reference_fasta, input_file = output_svsig, output_file = output_vcf)
		
		subprocess.run(pbsv_call_cmd, shell=True, check=True)

		compress_cmd = "bgzip -c {input_file} > {input_file}.gz".format(input_file = output_vcf)
		index_vcf_cmd = "tabix -p vcf {input_file}.gz".format(input_file = output_vcf)
		
		subprocess.run(compress_cmd, shell=True, check=True)
		subprocess.run(index_vcf_cmd, shell=True, check=True)

		print("SV VCF written to: {}".format(output_vcf))
		print("\n\n")

	# Genotype tandem repeats with pbtrgt
	def genotype_tandem_repeats(self):
		print("Genotyping tandem repeats with pbtrgt!")

		input_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")

		print("trgt input file: {}".format(input_bam))
		
		output_prefix = self.sample_ID + ".dedup.trimmed.hg38.chr6.TR"
		
		os.chdir(Samples.pbtrgt_dir)
		
		trgt_cmd = "trgt genotype --threads {threads} --genome {reference_genome} --reads {input_file} --repeats {repeat_file} --output-prefix {output_prefix} --preset targeted".format(threads = max_threads, reference_genome = reference_fasta, input_file = input_bam, repeat_file = pbtrgt_repeat_file, output_prefix = output_prefix)
		
		subprocess.run(trgt_cmd, shell=True, check=True)
		
		sort_cmd = "bcftools sort -O z -o {output_file} {input_file}".format(output_file = output_prefix + ".sorted.vcf.gz", input_file = output_prefix + ".vcf.gz")

		subprocess.run(sort_cmd, shell=True, check=True)

		os.rename(output_prefix + ".sorted.vcf.gz", output_prefix + ".vcf.gz")

		index_cmd = "tabix -p vcf {input_file}".format(input_file = output_prefix + ".vcf.gz")
		
		subprocess.run(index_cmd, shell=True, check=True)

		print("TR VCF written to {}".format(Samples.pbtrgt_dir + output_prefix + ".vcf.gz"))
		print("\n\n")

	# Merge SNV, tandem repeat, and structural variant vcfs with bcftools concat
	def merge_vcfs(self):
		print("Merging DeepVariant, pbsv, and pbtrgt VCF files!")

		input_snv = os.path.join(Samples.deepvariant_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.SNV.vcf.gz")
		input_SV = os.path.join(Samples.pbsv_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.SV.vcf.gz")
		input_TR = os.path.join(Samples.pbtrgt_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.TR.vcf.gz")

		print("DeepVariant input file: {}".format(input_snv))
		print("pbsv input file: {}".format(input_SV))
		print("pbtrgt input file: {}".format(input_TR))

		output_vcf = os.path.join(Samples.merged_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.vcf.gz")
		
		concat_cmd = "bcftools concat --allow-overlaps {SNV_vcf} {SV_vcf} {TR_vcf} | grep -v 'chrX|chrY' | grep -v 'SVTYPE=BND|SVTYPE=INV|SVTYPE=DUP' | bcftools norm -d none --fasta-ref {reference_genome} | bcftools sort | bgzip > {output_file}".format(SNV_vcf = input_snv, SV_vcf = input_SV, TR_vcf = input_TR, reference_genome = reference_fasta, output_file = output_vcf)
		
		subprocess.run(concat_cmd, shell=True, check=True)

		index_cmd = "tabix {output_file}".format(output_file = output_vcf)
		
		subprocess.run(index_cmd, shell=True, check=True)

		print("Merged VCF written to: {}".format(output_vcf))
		print("\n\n")

	# Phase genotypes with HiPhase
	def phase_genotypes_hiphase(self):
		print("Phasing Genotypes with HiPhase!")

		input_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")
		input_vcf = os.path.join(Samples.merged_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.vcf.gz")

		print("Input BAM: {}".format(input_bam))
		print("Input VCF: {}".format(input_vcf))
		
		output_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.haplotag.bam")
		output_vcf = os.path.join(Samples.hiphase_phased_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.phased.vcf.gz")
		output_summary_file = os.path.join(Samples.hiphase_phased_vcf_dir, self.sample_ID + ".phased.summary.txt")
		output_blocks_file = os.path.join(Samples.hiphase_phased_vcf_dir, self.sample_ID + ".phased.blocks.txt")
		output_stats_file = os.path.join(Samples.hiphase_phased_vcf_dir, self.sample_ID + ".phased.stats.txt")

		hiphase_cmd = "hiphase --threads {threads} --ignore-read-groups --reference {reference_genome} --bam {in_bam} --output-bam {out_bam} --vcf {in_vcf} --output-vcf {out_vcf} --stats-file {stats_file} --blocks-file {blocks_file} --summary-file {summary_file}".format(threads = max_threads, reference_genome = reference_fasta, in_bam = input_bam, out_bam = output_bam, in_vcf = input_vcf, out_vcf = output_vcf, stats_file = output_stats_file, blocks_file = output_blocks_file, summary_file = output_summary_file)
		
		# Log HiPhase in own output file so it doesn't clog up STDOUT
		hiphase_log = os.path.join(Samples.hiphase_phased_vcf_dir, self.sample_ID + ".hiphase.log")

		with open(hiphase_log, "w") as log_file:
			subprocess.run(hiphase_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

		print("HiPhase phased VCF written to: {}".format(output_vcf))
		print("HiPhase haplotagged BAM written to: {}".format(output_bam))
		print("HiPhase phasing summary written to: {}".format(output_summary_file))
		print("HiPhase phasing stats written to: {}".format(output_stats_file))
		print("HiPhase phase blocks written to: {}".format(output_blocks_file))
		print("\n\n")

	# Phase genotypes with WhatsHap
	def phase_genotypes_whatshap(self):
		print("Phasing Genotypes with WhatsHap!")

		input_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")
		input_vcf = os.path.join(Samples.deepvariant_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.SNV.vcf.gz")

		print("Input BAM: {}".format(input_bam))
		print("Input VCF: {}".format(input_vcf))

		haplotagged_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.haplotag.bam")
		phased_vcf = os.path.join(Samples.whatshap_phased_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.phased.vcf.gz")
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
		print("WhatsHap haplotagged BAM written to: {}".format(output_bam))
		print("WhatsHap phase block gtf written to: {}".format(output_gtf_file))
		print("WhatsHap phase blocks written to: {}".format(output_blocks_file))
		print("\n\n")

def main():
	# Check that all required tools are installed
	check_required_commands()

	for sample_ID, sample_read_group_string in sample_dict.items():
		start_time = time.time()
		sample = Samples(sample_ID, sample_read_group_string)
		# sample.convert_bam_to_fastq()
		# sample.mark_duplicates()
		# sample.run_fastqc(os.path.join(Samples.fastq_rmdup_dir, sample_ID + ".dedup.fastq.gz"))
		# sample.trim_adapters()
		# sample.run_fastqc(os.path.join(Samples.fastq_rmdup_cutadapt_dir, sample_ID + ".dedup.trimmed.fastq.gz"))
		# sample.align_to_reference()
		
		chr6_reads = sample.filter_reads()

		if chr6_reads > min_reads_sample:
			# sample.call_variants()
			# sample.call_structural_variants()
			# sample.genotype_tandem_repeats()
			# sample.merge_vcfs()
			# sample.phase_genotypes_hiphase()
			sample.phase_genotypes_whatshap()
			end_time = time.time()
			elapsed_time = end_time - start_time
			minutes, seconds = divmod(elapsed_time,60)
			print("Processed sampled in {}:{:.2f}!".format(int(minutes), seconds))
		
		else:
			print("Insufficient reads for variant calling")
			print("Sample {sample_id} had {num_reads} reads!".format(sample_id = sample_ID, num_reads = chr6_reads))

if __name__ == "__main__":
	main()
