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
mhc_classes = ["MHC_Class_I", "MHC_Class_II", "MHC_Class_III"]

root_dir = "/hb/scratch/mglasena/downsample_replicates/"
reference_fasta = os.path.join(root_dir, "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
deepvariant_sif = os.path.join(root_dir, "deepvariant_sif/deepvariant.sif")
temporary_directory = "/hb/scratch/mglasena/temp/"
mapped_bam_file = "/hb/groups/cornejo_lab/matt/pacbio_capture/processed_data/mapped_bam/HG002.dedup.trimmed.hg38.chr6.bam"
giab_benchmark_dir = "/hb/scratch/mglasena/MHC/concordance/GIAB_benchmark/"

regions_file_class_I = "/hb/scratch/mglasena/downsample_replicates/merged_hla_legacy_class_I_subset.bed"
regions_file_class_II = "/hb/scratch/mglasena/downsample_replicates/merged_hla_legacy_class_II_subset.bed"
regions_file_class_III = "/hb/scratch/mglasena/downsample_replicates/merged_hla_legacy_class_III_subset.bed"
rtg_path = "/hb/home/mglasena/.conda/envs/happy/bin/rtg"
rtg_template = "/hb/scratch/mglasena/MHC/concordance/hap_py_input/rtg_sdf_template"
tandem_repeat_bed = "/hb/scratch/mglasena/downsample_replicates/repeats_bed/human_GRCh38_no_alt_analysis_set.trf.bed"
pbtrgt_repeat_file = "/hb/scratch/mglasena/downsample_replicates/repeats_bed/polymorphic_repeats.hg38.bed"

mosdepth_regions_file_class_I = "/hb/scratch/mglasena/downsample_replicates/mhc_class_I_mosdepth_regions.bed"
mosdepth_regions_file_class_II = "/hb/scratch/mglasena/downsample_replicates/mhc_class_II_mosdepth_regions.bed"
mosdepth_regions_file_class_III = "/hb/scratch/mglasena/downsample_replicates/mhc_class_III_mosdepth_regions.bed"

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

num_replicates = 10

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

def downsample_bams(proportion, rep, mhc_class):
	replicate_dir = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"))
	os.makedirs(replicate_dir, exist_ok=True)

	input_file = os.path.join(root_dir, f"{mhc_class}.bam")
	output_bam = os.path.join(replicate_dir, f"{mhc_class}_downsampled.bam")
	seed = random_seed + rep

	downsample_cmd = (
		"gatk DownsampleSam -I {input_bam} -O {output_bam} "
		"-P {probability} -S HighAccuracy -R {random_seed} "
		"--CREATE_INDEX true --TMP_DIR {temp_dir} --VALIDATION_STRINGENCY LENIENT"
	).format(
		input_bam=input_file,
		output_bam=output_bam,
		probability=proportion,
		random_seed=seed,
		temp_dir=temporary_directory
	)

	subprocess.run(downsample_cmd, shell=True, check=True)

def merge_downsampled_bams(proportion, rep):
	# Define input and output paths
	input_bams = [
	os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), "MHC_Class_I_downsampled.bam"),
	os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), "MHC_Class_II_downsampled.bam"),
	os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), "MHC_Class_III_downsampled.bam"),
	]

	output_bam = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), "merged_downsampled.bam")

	input_file_string = ' '.join(input_bams)

	merge_cmd = "samtools merge -@ {threads} -o {output_file} {input_files}".format(threads = max_threads, output_file = output_bam, input_files = input_file_string)
	index_cmd = "samtools index {input_file}".format(input_file = output_bam)

	subprocess.run(merge_cmd, shell=True, check=True)
	subprocess.run(index_cmd, shell=True, check=True)

def run_mosdepth(proportion, mhc_class, rep):
	if mhc_class == "MHC_Class_I":
		regions = mosdepth_regions_file_class_I
	elif mhc_class == "MHC_Class_II":
		regions = mosdepth_regions_file_class_II
	elif mhc_class == "MHC_Class_III":
		regions = mosdepth_regions_file_class_III
	
	input_dir = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"))
	input_bam = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + "_downsampled.bam")
	os.chdir(input_dir)

	threads = 4 
		
	# --flag 3328 excludes duplicates and secondary/supplementary alignments
	mosdepth_cmd = "mosdepth --flag 3328 --mapq {mapq} --by {regions} -t {threads} {prefix} {input_bam}".format(mapq = mapq_threshold, regions = regions, threads = 4, prefix = mhc_class, input_bam = input_bam)
	subprocess.run(mosdepth_cmd, shell=True, check=True)

def parse_mosdepth(proportion, mhc_class, rep):
	replicate_id = f"{proportion}_rep{rep}"
	coverage_dict[mhc_class][replicate_id] = {}
	
	regions_file = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + ".regions.bed.gz")
	
	with gzip.open(regions_file, "rt") as f:
		regions = f.read().splitlines()

		for region in regions:
			fields = region.split("\t")
			gene = fields[3].split("_")[0]
			coverage_depth = float(fields[4])
			coverage_dict[mhc_class][replicate_id][gene] = coverage_depth

def count_reads(proportion, mhc_class, rep):
	#input_bam = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), "merged_downsampled.bam")
	input_bam = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + "_downsampled.bam")

	count_cmd = "samtools view -c {input_file}".format(input_file = input_bam)
	read_count = int(subprocess.check_output(count_cmd, shell=True).decode("utf-8").strip())

	return read_count

def call_variants(proportion, mhc_class, rep):
	input_dir = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"))
	#input_bam = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), "merged_downsampled.bam")
	input_bam = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + "_downsampled.bam")

	print("Running DeepVariant on {}".format(input_bam))

	index_cmd = f"samtools index {input_bam}"
	subprocess.run(index_cmd, shell=True, check=True)

	#output_vcf = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), "merged_downsampled.vcf.gz")
	#output_gvcf = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), "merged_downsampled.g.vcf.gz")
	output_vcf = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + "_downsampled.vcf.gz")
	output_gvcf = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + "_downsampled.g.vcf.gz")

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
	deepvariant_log = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), f"{mhc_class}_deepvariant.log")

	with open(deepvariant_log, "w") as log_file:
		subprocess.run(deepvariant_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	print("VCF written to {}".format(output_vcf))
	print("GVCF written to {}".format(output_gvcf))
	print("\n\n")

def run_happy(proportion, mhc_class, rep):
	outdir = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"))
	output_prefix = os.path.join(outdir, mhc_class)

	if mhc_class == "MHC_Class_I":
		regions = regions_file_class_I
	elif mhc_class == "MHC_Class_II":
		regions = regions_file_class_II
	elif mhc_class == "MHC_Class_III":
		regions = regions_file_class_III

	#query_vcf = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), "merged_downsampled.vcf.gz")
	query_vcf = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + "_downsampled.vcf.gz")

	print("Benchmarking {}".format(query_vcf))

	truth_vcf = os.path.join(giab_benchmark_dir, sample, f"{sample}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
	
	confident_regions = os.path.join(giab_benchmark_dir, sample, f"{sample}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed")

	run_happy = "hap.py {truth} {query} -f {confident_regions} -R {regions} -r {ref} -o {output_prefix} --engine vcfeval --engine-vcfeval-path {rtg_path} --engine-vcfeval-template {rtg_template} --threads {threads}".format(
		truth = truth_vcf, query = query_vcf, confident_regions = confident_regions, regions = regions, ref = reference_fasta, output_prefix = output_prefix, rtg_path = rtg_path, rtg_template = rtg_template, threads = max_threads)

	subprocess.run(run_happy, shell=True, check=True)

def parse_happy(proportion, mhc_class, rep):
	replicate_id = f"{proportion}_rep{rep}"

	happy_output = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), f"{mhc_class}.summary.csv")

	with open(happy_output, "r") as f:
		reader = csv.reader(f)
		headers = next(reader)

		# Get column indices
		type_idx = headers.index("Type")
		filter_idx = headers.index("Filter")
		recall_idx = headers.index("METRIC.Recall")
		precision_idx = headers.index("METRIC.Precision")
		f1_idx = headers.index("METRIC.F1_Score")

		concordance_dict[mhc_class][replicate_id] = {}

		# Read and process each row
		for row in reader:
			variant_type = row[type_idx]
			filter_value = row[filter_idx]
			
			# Only keep rows where Filter == "PASS"
			if filter_value == "PASS":
				recall = float(row[recall_idx])
				precision = float(row[precision_idx])
				f1_score = float(row[f1_idx])

				concordance_dict[mhc_class][replicate_id][variant_type] = [recall, precision, f1_score]

def write_results(mhc_class, rep):
	output_csv = f"{mhc_class}_rep{rep}_summary.csv"
	with open(output_csv, "w", newline="") as f:
		writer = csv.writer(f)
		header = ["Proportion", "Depth", "Variant", "Metric", "Value"]
		writer.writerow(header)

		for replicate_id, variants in concordance_dict[mhc_class].items():
			avg_depth = sum(coverage_dict[mhc_class][replicate_id].values()) / len(coverage_dict[mhc_class][replicate_id])

			for variant_type, metrics in variants.items():
				recall, precision, f1_score = metrics

				writer.writerow([replicate_id, avg_depth, variant_type, "Recall", recall])
				writer.writerow([replicate_id, avg_depth, variant_type, "Precision", precision])
				writer.writerow([replicate_id, avg_depth, variant_type, "F1", f1_score])

# Run pbsv to call structural variants (SV)
def call_structural_variants_pbsv(proportion, mhc_class, rep):
	print("Calling structural variants with pbsv!")

	input_bam = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + "_downsampled.bam")

	print("pbsv input file: {}".format(input_bam))

	# Ensure the BAM file is indexed as .bam.bai
	index_cmd = f"samtools index {input_bam}"
	subprocess.run(index_cmd, shell=True, check=True)
	
	output_svsig = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + ".dedup.trimmed.hg38.chr6.svsig.gz")
	output_vcf = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + ".dedup.trimmed.hg38.chr6.SV.vcf")

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

	print("pbsv SV VCF written to: {}".format(output_vcf))
	print("\n\n")
	# Genotype tandem repeats with pbtrgt
	
def genotype_tandem_repeats(proportion, mhc_class, rep):
	print("Genotyping tandem repeats with pbtrgt!")

	input_dir = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"))
	input_bam = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + "_downsampled.bam")

	print("trgt input file: {}".format(input_bam))
	
	output_prefix = mhc_class + ".dedup.trimmed.hg38.chr6.TR"
	
	os.chdir(input_dir)
	
	trgt_cmd = "trgt genotype --threads {threads} --genome {reference_genome} --reads {input_file} --repeats {repeat_file} --output-prefix {output_prefix} --preset targeted".format(threads = max_threads, reference_genome = reference_fasta, input_file = input_bam, repeat_file = pbtrgt_repeat_file, output_prefix = output_prefix)
	
	subprocess.run(trgt_cmd, shell=True, check=True)
	
	sort_cmd = "bcftools sort -O z -o {output_file} {input_file}".format(output_file = output_prefix + ".sorted.vcf.gz", input_file = output_prefix + ".vcf.gz")

	subprocess.run(sort_cmd, shell=True, check=True)

	os.rename(output_prefix + ".sorted.vcf.gz", output_prefix + ".vcf.gz")

	index_cmd = "tabix -p vcf {input_file}".format(input_file = output_prefix + ".vcf.gz")
	
	subprocess.run(index_cmd, shell=True, check=True)

	print("TR VCF written to {}".format(os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), output_prefix + ".vcf.gz")))
	print("\n\n")

# Merge SNV (DeepVariant), tandem repeat (pbtrgt), and structural variant (pbsv) vcfs with bcftools concat for use with HiPhase
def merge_vcfs(proportion, mhc_class, rep):
	print("Merging DeepVariant, pbsv, and pbtrgt VCF files!")

	input_snv = output_vcf = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + "_downsampled.vcf.gz")
	input_SV = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + ".dedup.trimmed.hg38.chr6.SV.vcf.gz")
	input_TR = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + ".dedup.trimmed.hg38.chr6.TR.vcf.gz")

	print("DeepVariant input file: {}".format(input_snv))
	print("pbsv input file: {}".format(input_SV))
	print("pbtrgt input file: {}".format(input_TR))

	output_vcf = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + ".dedup.trimmed.hg38.chr6.vcf.gz")
	
	concat_cmd = "bcftools concat --allow-overlaps {SNV_vcf} {SV_vcf} {TR_vcf} | grep -v -E 'chrX|chrY' | grep -v -E 'SVTYPE=BND|SVTYPE=INV|SVTYPE=DUP' | bcftools norm -d none --fasta-ref {reference_genome} | bcftools sort | bgzip > {output_file}".format(SNV_vcf = input_snv, SV_vcf = input_SV, TR_vcf = input_TR, reference_genome = reference_fasta, output_file = output_vcf)
	
	subprocess.run(concat_cmd, shell=True, check=True)

	index_cmd = "tabix {output_file}".format(output_file = output_vcf)
	
	subprocess.run(index_cmd, shell=True, check=True)

	print("Merged VCF written to: {}".format(output_vcf))
	print("\n\n")

# Phase genotypes with HiPhase
def phase_genotypes_hiphase(proportion, mhc_class, rep):
	print("Phasing Genotypes with HiPhase!")

	input_bam = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + "_downsampled.bam")
	input_vcf = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + ".dedup.trimmed.hg38.chr6.vcf.gz")

	print("Input BAM: {}".format(input_bam))
	print("Input VCF: {}".format(input_vcf))
	
	output_bam = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + ".dedup.trimmed.hg38.chr6.hiphase.haplotag.bam")
	output_vcf = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + ".dedup.trimmed.hg38.chr6.phased.vcf.gz")
	output_summary_file = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + ".phased.summary.txt")
	output_blocks_file = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + ".phased.blocks.txt")
	output_stats_file = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + ".phased.stats.txt")

	hiphase_cmd = "hiphase --threads {threads} --ignore-read-groups --reference {reference_genome} --bam {in_bam} --output-bam {out_bam} --vcf {in_vcf} --output-vcf {out_vcf} --stats-file {stats_file} --blocks-file {blocks_file} --summary-file {summary_file}".format(threads = max_threads, reference_genome = reference_fasta, in_bam = input_bam, out_bam = output_bam, in_vcf = input_vcf, out_vcf = output_vcf, stats_file = output_stats_file, blocks_file = output_blocks_file, summary_file = output_summary_file)
	
	# Log HiPhase in own output file so it doesn't clog up STDOUT
	hiphase_log = os.path.join(root_dir, os.path.join(str(proportion), f"rep{rep}"), mhc_class + ".hiphase.log")

	with open(hiphase_log, "w") as log_file:
		subprocess.run(hiphase_cmd, shell=True, check=True, stdout=log_file, stderr=log_file)

	print("HiPhase phased VCF written to: {}".format(output_vcf))
	print("HiPhase haplotagged BAM written to: {}".format(output_bam))
	print("HiPhase phasing summary written to: {}".format(output_summary_file))
	print("HiPhase phasing stats written to: {}".format(output_stats_file))
	print("HiPhase phase blocks written to: {}".format(output_blocks_file))
	print("\n\n")

def main():
	# make_output_dirs()
	# split_bam()
	array_id = int(os.environ["array_id"])
	print("Array ID: {}".format(array_id))
	proportion = proportion_retain[array_id]
	print(f"Processing proportion: {proportion}")

	for rep in range(num_replicates):
		# Downsample all three MHC class BAMs
		for mhc_class in mhc_classes:
			downsample_bams(proportion, rep, mhc_class)
		
		# Merge once all 3 are ready
		merge_downsampled_bams(proportion, rep)
		
		for mhc_class in mhc_classes:
			run_mosdepth(proportion, mhc_class, rep)
			parse_mosdepth(proportion, mhc_class, rep)

			num_reads = count_reads(proportion, mhc_class, rep)

			if num_reads >= min_reads:
				call_variants(proportion, mhc_class, rep)
				# run_happy(proportion, mhc_class, rep)
				# parse_happy(proportion, mhc_class, rep)
				# call_structural_variants_pbsv(proportion, mhc_class, rep)
				# genotype_tandem_repeats(proportion, mhc_class, rep)
				# merge_vcfs(proportion, mhc_class, rep)
				# phase_genotypes_hiphase(proportion, mhc_class, rep)
			else:
				print(f"Insufficient reads for variant calling at {proportion}_rep{rep} for {mhc_class}")

		# Write result *once per rep* after all classes processed
		for mhc_class in mhc_classes:
			write_results(mhc_class, rep)

if __name__ == "__main__":
	main()
