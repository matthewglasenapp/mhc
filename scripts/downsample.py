import subprocess
import os
import csv

# MHC Class I
MHC_I = "chr6:29555628-31511124"
# MHC Class III
MHC_III = "chr6:31519479-32407181"
# MHC Class II
MHC_II = "chr6:32439877-33409896"

sample = "HG002"

root_dir = "/hb/scratch/mglasena/downsample/"
reference_fasta = os.path.join(root_dir, "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
deepvariant_sif = os.path.join(root_dir, "deepvariant_sif/deepvariant.sif")
temporary_directory = "/hb/scratch/mglasena/temp/"
mapped_bam_file = "/hb/scratch/mglasena/test_pacbio/processed_data/mapped_bam/HG002.dedup.trimmed.hg38.chr6.bam"
giab_benchmark_dir = "/hb/scratch/mglasena/MHC/concordance/GIAB_benchmark/"
regions_file = "/hb/scratch/mglasena/MHC/concordance/hap_py_input/merged_hla_legacy.bed"
rtg_path = "/hb/home/mglasena/.conda/envs/happy/bin/rtg"
rtg_template = "/hb/scratch/mglasena/MHC/concordance/hap_py_input/rtg_sdf_template"

# Random Seed
random_seed = 42

# Threads
max_threads = 24

min_reads = 100

# Downsample proportions
proportion_retain = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
#proportion_retain = [0.5]

concordance_dict = dict()
output_csv = "downsample_concordance.csv"

def make_output_dirs():
	for n in proportion_retain:
		os.makedirs(os.path.join(root_dir, str(n)), exist_ok=True)

def split_bam():
	output_directory = root_dir
	subset_MHC_I = "samtools view -b {input_bam} {region} > {outdir}MHC_Class_I.bam".format(input_bam = mapped_bam_file, region = MHC_I, outdir = output_directory)
	subset_MHC_II = "samtools view -b {input_bam} {region} > {outdir}MHC_Class_II.bam".format(input_bam = mapped_bam_file, region = MHC_II, outdir = output_directory)
	subset_MHC_III = "samtools view -b {input_bam} {region} > {outdir}MHC_Class_III.bam".format(input_bam = mapped_bam_file, region = MHC_III, outdir = output_directory)

	subprocess.run(subset_MHC_I, shell=True, check=True)
	subprocess.run(subset_MHC_II, shell=True, check=True)
	subprocess.run(subset_MHC_III, shell=True, check=True)

def downsample_bams(proportion):
	input_BAMs = ["MHC_Class_I.bam", "MHC_Class_II.bam", "MHC_Class_III.bam"]

	for bam in input_BAMs:
		input_file = os.path.join(root_dir, bam)
		output_prefix = bam.split(".")[0] + "_downsampled.bam"
		output_bam = os.path.join(root_dir, str(proportion), output_prefix)
		downsample_cmd = "gatk DownsampleSam -I {input_bam} -O {output_bam} -P {probability} -S HighAccuracy -R {random_seed} --CREATE_INDEX true --TMP_DIR {temp_dir} --VALIDATION_STRINGENCY LENIENT".format(input_bam = input_file, output_bam = output_bam, probability = proportion, random_seed = random_seed, temp_dir = temporary_directory)

		subprocess.run(downsample_cmd, shell=True, check=True)

def merge_downsampled_bams(proportion):

	# Define input and output paths
	input_bams = [
	os.path.join(root_dir, str(proportion), "MHC_Class_I_downsampled.bam"),
	os.path.join(root_dir, str(proportion), "MHC_Class_II_downsampled.bam"),
	os.path.join(root_dir, str(proportion), "MHC_Class_III_downsampled.bam"),
	]

	output_bam = os.path.join(root_dir, str(proportion), "merged_downsampled.bam")

	input_file_string = ' '.join(input_bams)

	merge_cmd = "samtools merge -@ {threads} -o {output_file} {input_files}".format(threads = max_threads, output_file = output_bam, input_files = input_file_string)
	index_cmd = "samtools index {input_file}".format(input_file = output_bam)

	subprocess.run(merge_cmd, shell=True, check=True)
	subprocess.run(index_cmd, shell=True, check=True)

def count_reads(proportion):
	input_bam = os.path.join(root_dir, str(proportion), "merged_downsampled.bam")

	count_cmd = "samtools view -c {input_file}".format(input_file = input_bam)
	read_count = int(subprocess.check_output(count_cmd, shell=True).decode("utf-8").strip())

	return read_count

def call_variants(proportion):
	input_dir = os.path.join(root_dir, str(proportion))
	input_bam = os.path.join(root_dir, str(proportion), "merged_downsampled.bam")

	print("Running DeepVariant on {}".format(input_bam))

	output_vcf = os.path.join(root_dir, str(proportion), "merged_downsampled.vcf.gz")
	output_gvcf = os.path.join(root_dir, str(proportion), "merged_downsampled.g.vcf.gz")

	ref_filename = os.path.basename(reference_fasta)

	bind_paths = [
		f"{input_dir}:/data",
		f"{input_dir}:/input",
		f"{os.path.dirname(reference_fasta)}:/reference"
	]

	bind_flags = " ".join("--bind {}".format(path) for path in bind_paths)

	deepvariant_cmd = """
		singularity exec {binds} {sif} /opt/deepvariant/bin/run_deepvariant \
			--model_type=PACBIO \
			--ref=/reference/{ref_filename} \
			--reads=/input/merged_downsampled.bam \
			--output_vcf=/data/merged_downsampled.vcf.gz \
			--output_gvcf=/data/merged_downsampled.g.vcf.gz \
			--regions chr6 \
			--num_shards=8
		""".format(
			binds=bind_flags,
			sif=deepvariant_sif,
			ref_filename=os.path.basename(reference_fasta),
			sample=sample
			)

	# Log DeepVariant in own output file so it doesn't clog up STDOUT
	deepvariant_log = os.path.join(root_dir, str(proportion), "deepvariant.log")

	with open(deepvariant_log, "w") as log_file:
		subprocess.run(deepvariant_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	print("VCF written to {}".format(output_vcf))
	print("GVCF written to {}".format(output_gvcf))
	print("\n\n")

def run_happy(proportion):
	outdir = os.path.join(root_dir, str(proportion))
	output_prefix = os.path.join(outdir, sample)
	query_vcf = os.path.join(root_dir, str(proportion), "merged_downsampled.vcf.gz")

	print("Benchmarking {}".format(query_vcf))

	truth_vcf = os.path.join(giab_benchmark_dir, sample, f"{sample}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
	
	confident_regions = os.path.join(giab_benchmark_dir, sample, f"{sample}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed")

	run_happy = "hap.py {truth} {query} -f {confident_regions} -R {regions} -r {ref} -o {output_prefix} --engine vcfeval --engine-vcfeval-path {rtg_path} --engine-vcfeval-template {rtg_template} --threads {threads}".format(
		truth = truth_vcf, query = query_vcf, confident_regions = confident_regions, regions = regions_file, ref = reference_fasta, output_prefix = output_prefix, rtg_path = rtg_path, rtg_template = rtg_template, threads = max_threads)

	subprocess.run(run_happy, shell=True, check=True)

def parse_happy(proportion):
	happy_output = os.path.join(root_dir, str(proportion), f"{sample}.summary.csv")

	with open(happy_output, "r") as f:
		reader = csv.reader(f)
		headers = next(reader)

		# Get column indices
		type_idx = headers.index("Type")
		filter_idx = headers.index("Filter")
		recall_idx = headers.index("METRIC.Recall")
		precision_idx = headers.index("METRIC.Precision")
		f1_idx = headers.index("METRIC.F1_Score")

		concordance_dict[str(proportion)] = {}

		# Read and process each row
		for row in reader:
			variant_type = row[type_idx]
			filter_value = row[filter_idx]
			
			# Only keep rows where Filter == "PASS"
			if filter_value == "PASS":
				recall = float(row[recall_idx])
				precision = float(row[precision_idx])
				f1_score = float(row[f1_idx])

				concordance_dict[str(proportion)][variant_type] = [recall, precision, f1_score]

def write_results():
	with open(output_csv, "w", newline="") as f:
		writer = csv.writer(f)
		header = ["Proportion", "Variant", "Metric", "Value"]
		writer.writerow(header)

		for proportion, variants in concordance_dict.items():
			for variant_type, metrics in variants.items():
				recall, precision, f1_score = metrics

				writer.writerow([proportion, variant_type, "Recall", recall])
				writer.writerow([proportion, variant_type, "Precision", precision])
				writer.writerow([proportion, variant_type, "F1", f1_score])

def main():
	make_output_dirs()
	split_bam()
	for n in proportion_retain:
		downsample_bams(n)
		merge_downsampled_bams(n)

		num_reads = count_reads(n)

		if num_reads >= min_reads:
			call_variants(n)
			#run_happy(n)
			#parse_happy(n)
		else:
			print("Insufficient Reads for Variant Calling at {} Downample".format(n))

	#write_results()

if __name__ == "__main__":
	main()
