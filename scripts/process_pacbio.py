import os
import shutil
import subprocess
import sys

"""
Work Flow

1. bam2fastq
2. pbmarkdup
3. cutadapt
4. pbmm2
5. DeepVariant
6. pbsv discover
7. pbsv call
8. pbtrgt
9. HiPhase
"""

# Try to run with as few compute resources as possible and document
# Document the size of data (average per sample) at each step.
# One sample takes X amount of time on Y core with Z CPU/RAM 
# Run time python3 -u process_pacbio.py with 24 CPU, 60 GB RAM

# set input directory to the current working directory where the script should be run
input_dir = os.getcwd()

# Set the output directory to processed_data/
output_dir = os.path.join(input_dir, "processed_data")
os.makedirs(output_dir, exist_ok=True)

# Input file paths
# Someone reccomends using fasta with no alternate contigs.
# GRCh38 tandem repeat bed file for pbsv
# Downloaded from https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
# Repeat definition file for pbtrgt
# Downloaded from https://zenodo.org/records/8329210
reference_fasta = os.path.join(input_dir, "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
deepvariant_sif = os.path.join(input_dir, "deepvariant_sif/deepvariant.sif")
tandem_repeat_bed = os.path.join(input_dir, "repeats_bed/human_GRCh38_no_alt_analysis_set.trf.bed")
pbtrgt_repeat_file = os.path.join(input_dir, "repeats_bed/polymorphic_repeats.hg38.bed")

# Sample ID: [Read Group String, Karyotype]
sample_dict = {
"HG002": [r"@RG\tID:m84039_240622_113450_s1\tSM:HG002", "XY"]
}

# Max threads available for parallelization
max_threads = 8

# Transposase mosaic end binding sequence
# Adapters and barcodes were removed by PacBio with lima
# The TE sequence (and its reverse complement) introduced during tagmentation still needs to be removed
me = "AGATGTGTATAAGAGACAG"
me_rc = "CTGTCTCTTATACACATCT"

def check_required_commands():    
	print("Checking the installation status of the required bioinformatics tools!")

	required_commands = [
		"bam2fastq",
		"pbmarkdup",
		"pigz",
		"fastqc",
		"cutadapt",
		"pbmm2",
		"samtools",
		"singularity",
		"pbsv",
		"trgt",
		"bcftools",
		"hiphase"
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
	raw_hifi_reads_dir = os.path.join(output_dir, "raw_hifi_reads")
	fastq_raw_dir = os.path.join(output_dir, "fastq_raw")
	fastq_rmdup_dir = os.path.join(output_dir, "fastq_rmdup")
	fastq_rmdup_cutadapt_dir = os.path.join(output_dir, "fastq_rmdup_cutadapt")
	mapped_bam_dir = os.path.join(output_dir, "mapped_bam")
	deepvariant_dir = os.path.join(output_dir, "deepvariant_vcf")
	pbsv_dir = os.path.join(output_dir, "pbsv_vcf")
	pbtrgt_dir = os.path.join(output_dir, "pbtrgt_vcf")
	merged_vcf_dir = os.path.join(output_dir, "merged_vcf")
	phased_vcf_dir = os.path.join(output_dir, "phased_vcf")

	def __init__(self, sample_ID, read_group_string, karyotype):
		self.sample_ID = sample_ID
		self.unmapped_bam = os.path.join(input_dir, "raw_hifi_reads", self.sample_ID + ".hifi_reads.bam")
		self.read_group_string = read_group_string
		self.sample_karyotype = karyotype

		for directory in [Samples.fastq_raw_dir, Samples.fastq_rmdup_dir, Samples.fastq_rmdup_cutadapt_dir, Samples.mapped_bam_dir, Samples.deepvariant_dir, Samples.pbsv_dir, Samples.pbtrgt_dir, Samples.merged_vcf_dir, Samples.phased_vcf_dir]:
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
		
		print("Filtered BAM records written to: {}".format(output_bam))
		print("\n\n")

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
		
		trgt_cmd = "trgt genotype --threads {threads} --genome {reference_genome} --reads {input_file} --repeats {repeat_file} --output-prefix {output_prefix} --karyotype {karyotype} --preset targeted".format(threads = max_threads, reference_genome = reference_fasta, input_file = input_bam, repeat_file = pbtrgt_repeat_file, output_prefix = output_prefix, karyotype = self.sample_karyotype)
		
		#subprocess.run(trgt_cmd, shell=True, check=True)
		
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
	def phase_genotypes(self):
		print("Phasing Genotypes with HiPhase!")

		input_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.bam")
		input_vcf = os.path.join(Samples.merged_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.vcf.gz")

		print("Input BAM: {}".format(input_bam))
		print("Input VCF: {}".format(input_vcf))
		
		output_bam = os.path.join(Samples.mapped_bam_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.haplotag.bam")
		output_vcf = os.path.join(Samples.phased_vcf_dir, self.sample_ID + ".dedup.trimmed.hg38.chr6.phased.vcf.gz")
		output_summary_file = os.path.join(Samples.phased_vcf_dir, self.sample_ID + ".phased.summary.txt")
		output_blocks_file = os.path.join(Samples.phased_vcf_dir, self.sample_ID + ".phased.blocks.txt")
		output_stats_file = os.path.join(Samples.phased_vcf_dir, self.sample_ID + ".phased.stats.txt")

		hiphase_cmd = "hiphase --threads {threads} --ignore-read-groups --reference {reference_genome} --bam {in_bam} --output-bam {out_bam} --vcf {in_vcf} --output-vcf {out_vcf} --stats-file {stats_file} --blocks-file {blocks_file} --summary-file {summary_file}".format(threads = max_threads, reference_genome = reference_fasta, in_bam = input_bam, out_bam = output_bam, in_vcf = input_vcf, out_vcf = output_vcf, stats_file = output_stats_file, blocks_file = output_blocks_file, summary_file = output_summary_file)
		
		hiphase_log = os.path.join(Samples.phased_vcf_dir, self.sample_ID + ".hiphase.log")

		with open(hiphase_log, "w") as log_file:
			subprocess.run(hiphase_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

		print("Phased VCF written to: {}".format(output_vcf))
		print("Haplotagged BAM written to: {}".format(output_bam))
		print("Phasing summary written to: {}".format(output_summary_file))
		print("Phasing stats written to: {}".format(output_stats_file))
		print("Phase blocks written to: {}".format(output_blocks_file))
		print("\n\n")

def main():
	# Check that all required tools are installed
	check_required_commands()

	sample_ID = "HG002"
	sample = Samples(sample_ID, sample_dict[sample_ID][0], sample_dict[sample_ID][1])
	sample.convert_bam_to_fastq()
	sample.mark_duplicates()
	sample.run_fastqc(os.path.join(Samples.fastq_rmdup_dir, sample_ID + ".dedup.fastq.gz"))
	sample.trim_adapters()
	sample.run_fastqc(os.path.join(Samples.fastq_rmdup_cutadapt_dir, sample_ID + ".dedup.trimmed.fastq.gz"))
	sample.align_to_reference()
	sample.filter_reads()
	sample.call_variants()
	sample.call_structural_variants()
	sample.genotype_tandem_repeats()
	sample.merge_vcfs()
	sample.phase_genotypes()

if __name__ == "__main__":
	main()
