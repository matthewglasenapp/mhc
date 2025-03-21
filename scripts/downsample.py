import subprocess
import os
import csv
import gzip

# MHC Class I
MHC_I = "chr6:29555628-31511124"
# MHC Class III
MHC_III = "chr6:31519479-32407181"
# MHC Class II
MHC_II = "chr6:32439877-33409896"

sample = "HG002"
mhc_classes = ["MHC_Class_I", "MHC_Class_II"]

root_dir = "/hb/scratch/mglasena/downsample/"
reference_fasta = os.path.join(root_dir, "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
deepvariant_sif = os.path.join(root_dir, "deepvariant_sif/deepvariant.sif")
temporary_directory = "/hb/scratch/mglasena/temp/"
mapped_bam_file = "/hb/scratch/mglasena/test_pacbio/processed_data/mapped_bam/HG002.dedup.trimmed.hg38.chr6.bam"
giab_benchmark_dir = "/hb/scratch/mglasena/MHC/concordance/GIAB_benchmark/"
regions_file_class_I = "/hb/scratch/mglasena/downsample/merged_hla_legacy_class_I_subset.bed"
regions_file_class_II = "/hb/scratch/mglasena/downsample/merged_hla_legacy_class_II_subset.bed"
rtg_path = "/hb/home/mglasena/.conda/envs/happy/bin/rtg"
rtg_template = "/hb/scratch/mglasena/MHC/concordance/hap_py_input/rtg_sdf_template"

mosdepth_regions_file_class_I = "/hb/scratch/mglasena/downsample/mhc_class_I_mosdepth_regions.bed"
mosdepth_regions_file_class_II = "/hb/scratch/mglasena/downsample/mhc_class_II_mosdepth_regions.bed"

# Random Seed
random_seed = 42

# Threads
max_threads = 24

# Read threshold for DeepVariant
min_reads = 10

# MapQ threshold for mosdepth
mapq_threshold = 20

# Downsample proportions
proportion_retain = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 
 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1]

concordance_dict = {mhc_class: {} for mhc_class in mhc_classes}

coverage_dict = {mhc_class: {} for mhc_class in mhc_classes}

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

def run_mosdepth(proportion, mhc_class):
	if mhc_class == "MHC_Class_I":
		regions = mosdepth_regions_file_class_I
	elif mhc_class == "MHC_Class_II":
		regions = mosdepth_regions_file_class_II
	
	input_dir = os.path.join(root_dir, str(proportion))
	input_bam = os.path.join(root_dir, str(proportion), mhc_class + "_downsampled.bam")
	os.chdir(input_dir)

	threads = 4 
		
	# --flag 3328 excludes duplicates and secondary/supplementary alignments
	mosdepth_cmd = "mosdepth --flag 3328 --mapq {mapq} --by {regions} -t {threads} {prefix} {input_bam}".format(mapq = mapq_threshold, regions = regions, threads = 4, prefix = mhc_class, input_bam = input_bam)
	subprocess.run(mosdepth_cmd, shell=True, check=True)

def parse_mosdepth(proportion, mhc_class):
	coverage_dict[mhc_class][str(proportion)] = {}
	
	regions_file = os.path.join(root_dir, str(proportion), mhc_class + ".regions.bed.gz")
	
	with gzip.open(regions_file, "rt") as f:
		regions = f.read().splitlines()

		for region in regions:
			fields = region.split("\t")
			gene = fields[3].split("_")[0]
			coverage_depth = float(fields[4])
			coverage_dict[mhc_class][str(proportion)][gene] = coverage_depth

def count_reads(proportion, mhc_class):
	#input_bam = os.path.join(root_dir, str(proportion), "merged_downsampled.bam")
	input_bam = os.path.join(root_dir, str(proportion), mhc_class + "_downsampled.bam")

	count_cmd = "samtools view -c {input_file}".format(input_file = input_bam)
	read_count = int(subprocess.check_output(count_cmd, shell=True).decode("utf-8").strip())

	return read_count

def call_variants(proportion, mhc_class):
	input_dir = os.path.join(root_dir, str(proportion))
	#input_bam = os.path.join(root_dir, str(proportion), "merged_downsampled.bam")
	input_bam = os.path.join(root_dir, str(proportion), mhc_class + "_downsampled.bam")

	print("Running DeepVariant on {}".format(input_bam))

	#output_vcf = os.path.join(root_dir, str(proportion), "merged_downsampled.vcf.gz")
	#output_gvcf = os.path.join(root_dir, str(proportion), "merged_downsampled.g.vcf.gz")
	output_vcf = os.path.join(root_dir, str(proportion), mhc_class + "_downsampled.vcf.gz")
	output_gvcf = os.path.join(root_dir, str(proportion), mhc_class + "_downsampled.g.vcf.gz")

	ref_filename = os.path.basename(reference_fasta)

	bind_paths = [
		f"{input_dir}:/data",
		f"{input_dir}:/input",
		f"{os.path.dirname(reference_fasta)}:/reference"
	]

	bind_flags = " ".join("--bind {}".format(path) for path in bind_paths)

	# deepvariant_cmd = """
	# 	singularity exec {binds} {sif} /opt/deepvariant/bin/run_deepvariant \
	# 		--model_type=PACBIO \
	# 		--ref=/reference/{ref_filename} \
	# 		--reads=/input/merged_downsampled.bam \
	# 		--output_vcf=/data/merged_downsampled.vcf.gz \
	# 		--output_gvcf=/data/merged_downsampled.g.vcf.gz \
	# 		--regions chr6 \
	# 		--num_shards=8
	# 	""".format(
	# 		binds=bind_flags,
	# 		sif=deepvariant_sif,
	# 		ref_filename=os.path.basename(reference_fasta),
	# 		sample=sample
	# 		)
	
	deepvariant_cmd = """
		singularity exec {binds} {sif} /opt/deepvariant/bin/run_deepvariant \
			--model_type=PACBIO \
			--ref=/reference/{ref_filename} \
			--reads=/input/{input_bam} \
			--output_vcf=/data/{output_vcf} \
			--output_gvcf=/data/{output_gvcf} \
			--regions chr6 \
			--num_shards=8
	""".format(
		binds=bind_flags,
		sif=deepvariant_sif,
		ref_filename=ref_filename,
		input_bam=os.path.basename(input_bam),
		output_vcf=os.path.basename(output_vcf),
		output_gvcf=os.path.basename(output_gvcf)
	)


	# Log DeepVariant in own output file so it doesn't clog up STDOUT
	deepvariant_log = os.path.join(root_dir, str(proportion), f"{mhc_class}_deepvariant.log")

	with open(deepvariant_log, "w") as log_file:
		subprocess.run(deepvariant_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	print("VCF written to {}".format(output_vcf))
	print("GVCF written to {}".format(output_gvcf))
	print("\n\n")

def run_happy(proportion, mhc_class):
	outdir = os.path.join(root_dir, str(proportion))
	output_prefix = os.path.join(outdir, mhc_class)

	if mhc_class == "MHC_Class_I":
		regions = regions_file_class_I
	elif mhc_class == "MHC_Class_II":
		regions = regions_file_class_II

	#query_vcf = os.path.join(root_dir, str(proportion), "merged_downsampled.vcf.gz")
	query_vcf = os.path.join(root_dir, str(proportion), mhc_class + "_downsampled.vcf.gz")

	print("Benchmarking {}".format(query_vcf))

	truth_vcf = os.path.join(giab_benchmark_dir, sample, f"{sample}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
	
	confident_regions = os.path.join(giab_benchmark_dir, sample, f"{sample}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed")

	run_happy = "hap.py {truth} {query} -f {confident_regions} -R {regions} -r {ref} -o {output_prefix} --engine vcfeval --engine-vcfeval-path {rtg_path} --engine-vcfeval-template {rtg_template} --threads {threads}".format(
		truth = truth_vcf, query = query_vcf, confident_regions = confident_regions, regions = regions, ref = reference_fasta, output_prefix = output_prefix, rtg_path = rtg_path, rtg_template = rtg_template, threads = max_threads)

	subprocess.run(run_happy, shell=True, check=True)

def parse_happy(proportion, mhc_class):
	happy_output = os.path.join(root_dir, str(proportion), f"{mhc_class}.summary.csv")

	with open(happy_output, "r") as f:
		reader = csv.reader(f)
		headers = next(reader)

		# Get column indices
		type_idx = headers.index("Type")
		filter_idx = headers.index("Filter")
		recall_idx = headers.index("METRIC.Recall")
		precision_idx = headers.index("METRIC.Precision")
		f1_idx = headers.index("METRIC.F1_Score")

		concordance_dict[mhc_class][str(proportion)] = {}

		# Read and process each row
		for row in reader:
			variant_type = row[type_idx]
			filter_value = row[filter_idx]
			
			# Only keep rows where Filter == "PASS"
			if filter_value == "PASS":
				recall = float(row[recall_idx])
				precision = float(row[precision_idx])
				f1_score = float(row[f1_idx])

				concordance_dict[mhc_class][str(proportion)][variant_type] = [recall, precision, f1_score]

def write_results(mhc_class):
	output_csv = mhc_class + "_downsample_concordance.csv"
	with open(output_csv, "w", newline="") as f:
		writer = csv.writer(f)
		header = ["Proportion", "Depth", "Variant", "Metric", "Value"]
		writer.writerow(header)

		for proportion, variants in concordance_dict[mhc_class].items():
			avg_depth = sum(coverage_dict[mhc_class][proportion].values()) / len(coverage_dict[mhc_class][proportion])

			for variant_type, metrics in variants.items():
				recall, precision, f1_score = metrics

				writer.writerow([proportion, avg_depth, variant_type, "Recall", recall])
				writer.writerow([proportion, avg_depth, variant_type, "Precision", precision])
				writer.writerow([proportion, avg_depth, variant_type, "F1", f1_score])

def main():
	make_output_dirs()
	split_bam()
	for mhc_class in mhc_classes:
		for n in proportion_retain:
			print(f"\nProcessing {mhc_class} at {n*100:.1f}% downsampling")
			downsample_bams(n)
			#merge_downsampled_bams(n)
			run_mosdepth(n, mhc_class)
			parse_mosdepth(n, mhc_class)

			num_reads = count_reads(n, mhc_class)

			if num_reads >= min_reads:
				call_variants(n, mhc_class)
				#run_happy(n, mhc_class)
				#parse_happy(n, mhc_class)
			else:
				print("Insufficient Reads for Variant Calling at {} Downample".format(n))

			#write_results(mhc_class)

if __name__ == "__main__":
	main()
