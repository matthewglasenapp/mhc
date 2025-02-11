import os

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

# To Do 
# After pbmm2, use samtools to throw away everything not aligned to chromosome 6.  
# Document the size of data (average per sample) at each step.
# One sample takes X amount of time on Y core with Z CPU/RAM 

# Max threads available for parallelization
max_threads = 24

# Sample ID: [Read Group String, Karyotype]
sample_dict = {
"HG002": ["@RG\tID:m84039_240622_113450_s1\tSM:HG002", "XY"]
}

working_directory = "/hb/scratch/mglasena/test_pacbio/"

# Someone reccomends using fasta with no alternate contigs. 
reference_fasta = working_directory + "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

deepvariant_sif = working_directory + "deepvariant_sif/deepvariant.sif"

# GRCh38 tandem repeat bed file for pbsv
# Downloaded from https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
tandem_repeat_bed = working_directory + "input_files/human_GRCh38_no_alt_analysis_set.trf.bed"

# Repeat definition file for pbtrgt
# Downloaded from https://zenodo.org/records/8329210
pbtrgt_repeat_file = working_directory + "input_files/polymorphic_repeats.hg38.bed"

# Transposase mosaic end binding sequence
# Adapters and barcodes were removed by PacBio with lima
# The TE sequence (and its reverse complement) introduced during tagmentation still needs to be removed
me = "AGATGTGTATAAGAGACAG"
me_rc = "CTGTCTCTTATACACATCT"

class Samples:

    def __init__(self, sample_ID, read_group_string, karyotype):
        self.sample_ID = sample_ID
        self.unampped_bam = working_directory = "unmapped_bam/" + "sample.hifi_reads.bam"
        self.read_group_string = read_group_string
        self.sample_karyotype = karyotype

	# Convert BAM file of unampped HiFi (ccs) reads to FASTQ format for marking duplicates and trimming adapters
	def convert_bam_to_fastq(self):
		bam2fastq_cmd = "bam2fastq -j {threads} {input_file} -o {output_prefix}".format(threads = 24, input_file = self.unmapped_bam, output_prefix = self.sample_ID)
		os.system(bam2fastq_cmd)

	# Mark PCR duplicates with pbmarkdup
	def mark_duplicates(self):
		input_fastq = self.sample_ID + ".fastq.gz"
		output_fastq = self.sample_ID + ".dedup.fastq.gz"
		pbmarkdup_cmd = "pbmarkdup -j {max_threads} --rmdup {input_file} {output_file}".format(threads = max_threads, input_file = input_fastq, output_file = output_fastq)
		os.system(pbmarkdup_cmd)

	# Run fastqc on fastq file
	def run_fastqc(self, fastqc_file):
		fastqc_cmd = "fastqc --memory {memory} {input_fastq}".format(memory = 5000, input_fastq = fastqc_file)
		os.system(fastqc_cmd)

	# Trim Adapter Sequences and polyA tails
	# Customize this command based on the fastqc output
	def trim_adapters(self):
		input_file = self.sample_ID + ".dedup.fastq.gz"
		output_file = self.sample_ID + ".dedup.trimmed.fastq.gz"
		
		# The following argument is used to trim polyA tails: -a 'A{{10}}N{{90}}. The double brackets are for python string interpolation
		# me and me_rc are the transposase mosaic end binding sequences defined at the top of the script
		cutadapt_cmd = "cutadapt -j {threads} -n {num_cuts_allowed} -g {five_prime_adapter} -a {three_prime_adapter} -a 'A{{10}}N{{90}} -o {output_fastq} {input_fastq}".format(threads = max_threads, num_cuts_allowed = 3, five_prime_adapter = me, three_prime_adapter = me_rc, output_fastq = output_file, input_fastq = input_file)

		os.system(cutadapt_cmd)

	# Align to GRCh38 reference genome with pbmm2
	def align_to_reference(self)
		input_fastq = self.sample_ID + ".dedup.trimmed.fastq.gz"
		output_bam = self.sample_ID + ".dedup.trimed.hg38.bam"

		pbmm2_cmd = "pbmm2 align -j {threads} {reference_genome} {input_file} {output_file} --sort --log-level INFO --unmapped --bam-index BAI --rg {rg_string}".format(threads = max_threads, reference_genome = reference_fasta, input_file = input_fastq, output_file = output_bam, rg_string = self.read_group_string)

	# Call SNV with DeepVariant
	def call_variants(self):
	    input_bam = "{}.dedup.trimmed.hg38.bam".format(self.sample_ID)
	    output_vcf = "{}.dedup.trimmed.hg38.SNV.vcf.gz".format(self.sample_ID)
	    output_gvcf = "{}.dedup.trimmed.hg38.SNV.g.vcf.gz".format(self.sample_ID)

	    deepvariant_cmd = "singularity exec --bind {working_dir}:/data {sif} /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref=/data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa --reads=/data/{input_bam} --output_vcf=/data/{output_vcf} --output_gvcf=/data/{output_gvcf} --num_shards=8".format(
	        working_dir=working_dir, sif=deepvariant_sif, input_bam=input_bam, output_vcf=output_vcf, output_gvcf=output_gvcf)

	    os.system(deepvariant_cmd)

	# Run pbsv to call structural variants (SV)
	def call_structural_variants(self):
		input_bam = self.sample_ID + ".dedup.trimmed.hg38.bam"
		output_svsig = self.sample_ID + ".dedup.trimmed.hg38.svsig.gz"
		output_vcf = self.sample_ID + ".dedup.trimmed.hg38.SV.vcf"

		pbsv_discover_cmd = "pbsv discover --region chr6 --tandem-repeats {tandem_repeat_file} {input_file} {output_file}".format(tandem_repeat_file = tandem_repeat_bed, input_file = input_bam, output_file = output_svsig)
		os.system(pbsv_discover_cmd)

		index_svsig_cmd = "tabix -c '#' -s 3 -b 4 -e 4 {svsig_file}".format(svsig_file = output_svsig)
		os.system(index_svsig_cmd)

		pbsv_call_cmd = "pbsv call -j {threads} --region chr6 --hifi {reference_genome} {input_file} {output_file}".format(threads = max_threads, input_file = output_svsig, output_file = output_vcf)

		compress_cmd = "bgzip -c {input_file}".format(input_file = output_vcf)
		index_vcf_cmd = "tabix -p vcf {input_file}".format(input_file = output_vcf + ".gz")

	# Genotype tandem repeats with pbtrgt
	def genotype_tandem_repeats(self):
		input_bam = self.sample_ID + ".dedup.trimmed.hg38.bam"
		output_prefix = self.sample_ID + ".dedup.trimmed.hg38.TR"
		trgt_cmd = "trgt genotype --threads {threads} --genome {reference_genome} --reads {input_file} --repeats {repeat_file} --output-prefix {output_prefix} --karyotype {karyotype} --preset targeted".format(threads = max_threads, reference_genome = reference_fasta, input_file = bam, output_prefix = self.sample_ID + ".trgt", karyotype = self.sample_karyotype)
		os.system(trgt_cmd)

	# Merge SNV, tandem repeat, and structural variant vcfs with bcftools concat
	def merge_vcfs(self):
		input_snv = self.sample_ID + ".dedup.trimmed.hg38.SNV.vcf.gz"
		input_SV = self.sample_ID + ".dedup.trimmed.hg38.SV.vcf.gz"
		input_TR = self.sample_ID + ".dedup.trimmed.hg38.TR.vcf.gz"
		output_vcf = self.sample_ID ".dedup.trimmed.hg38.vcf.gz"
		concat_cmd = "bcftools concat --allow-overlaps {SNV_vcf} {SV_vcf} {TR_vcf} | grep -v 'chrX|chrY' | grep -v 'SVTYPE=BND|SVTYPE=INV|SVTYPE=DUP' | bcftools norm -D --fasta-ref {reference_genome} | bcftools sort | bgzip > {output_file}".format(SNV_vcf = input_snv, SV_vcf = input_SV, TR_vcf = input_TR, reference_genome = reference_fasta, output_file = output_vcf)
		os.system(merge_vcf)

		index_cmd = "tabix {output_file}".format(output_file = output_vcf)
		os.system(index_cmd)

	# Phase genotypes with HiPhase
	def phase_genotypes(self):
		input_bam = self.sample_ID + ".dedup.trimmed.hg38.bam"
		input_vcf = self.sample_ID + ".dedup.trimmed.hg38.vcf.gz"
		output_bam = self.sample_ID + ".dedup.trimmed.hg38.haplotag.bam"
		output_vcf = self.sample_ID + ".dedup.trimmed.hg38.phased.vcf.gz"
		hiphase_cmd = "hiphase --threads {threads} --ignore-read-groups --reference {reference_genome} --bam {in_bam} --output-bam {out_bam} --vcf {in_vcf} --output-vcf {out_vcf} --stats-file {stats_file} --blocks-file {blocks_file} --summary-file {summary_file}".format(threads = max_threads, reference_genome = reference_fasta, in_bam = input_bam, out_bam = output_bam, in_vcf = input_vcf, out_vcf = output_vcf, stats_file = self.sample_ID ".phased.stats.txt", blocks_file = self.sample_ID + ".phased.blocks.txt", summary_file = self.sample_ID ".phased.summary.txt")
		os.system(hiphase_cmd)

def main():
	sample_ID = "HG002"
	sample = Samples(sample_ID, sample_dict[sample_ID][0], sample_dict[sample_ID][1])
	sample.convert_bam_to_fastq()
	sample.mark_duplicates()
	sample.run_fastqc(sample_ID + ".dedup.fastq.gz")
	sample.trim_adapters(sample_ID + ".dedup.trimmed.fastq.gz")
	sample.run_fastqc(fastq_file)
	sample.align_to_reference()
	sample.call_variants()
	sample.call_structural_variants()
	sample.genotype_tandem_repeats()
	sample.merge_vcfs()
	sample.phase_genotypes()


if __name__ == "__main__":
	main()




