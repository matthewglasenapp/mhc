import os
import csv

threads = 8

samples = ["HG002", "HG003", "HG004", "HG005"]
platforms = ["revio", "promethion"]
types = ["snp", "indel"]

root_dir = "/hb/scratch/mglasena/MHC/concordance/"
input_vcf_dir = "/hb/scratch/mglasena/MHC/genotypes/"
hap_py_input_dir = root_dir + "hap_py_input/"
giab_benchmark_dir = root_dir + "GIAB_benchmark/"
output_dir = root_dir + "hap_py_results/"
reference_fasta = hap_py_input_dir + "Homo_sapiens.GRCh38.dna.primary_assembly_renamed.fa"
#regions_file = hap_py_input_dir + "merged_hla_legacy.bed"
regions_file = hap_py_input_dir + "c4b.bed"

# Path to rtg tools
rtg_path = "/hb/home/mglasena/.conda/envs/happy/bin/rtg"
rtg_template = hap_py_input_dir + "rtg_sdf_template"

# Set HGREF environment variable (hap.py expects this)
os.environ['HGREF'] = reference_fasta

# Create output directory if it doesn't exist
make_output_dir = "mkdir -p {}".format(output_dir)
os.system(make_output_dir)

concordance_dict = {sample: {platform: {type_: [] for type_ in types} for platform in platforms} for sample in samples}
slurm_output_file = "concordance.out"
parsed_results_file = output_dir + "concordance_results.csv"

def run_happy(platform, sample):
	print("Benchmarking {} ({})!".format(sample, platform))

	outdir = output_dir + platform + "/"
	os.system("mkdir -p {}".format(outdir))
	output_prefix = outdir + platform + "_" + sample
	if platform == "revio":
		query_vcf = input_vcf_dir + platform + "/" + sample + ".hg38." + platform + ".vcf.gz"
	elif platform == "promethion":
		query_vcf = input_vcf_dir + platform + "/" + sample + "/merge_output.vcf.gz"

	truth_vcf = giab_benchmark_dir + sample + "/" + sample + "_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

	if sample == "HG002" or sample == "HG003" or sample == "HG004":
		confident_regions = giab_benchmark_dir + sample + "/" + sample + "_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
	else:
		confident_regions = giab_benchmark_dir + sample + "/" + sample + "_GRCh38_1_22_v4.2.1_benchmark.bed"

	run_happy = "hap.py {} {} -f {} -R {} -r {} -o {} --engine vcfeval --engine-vcfeval-path {} --engine-vcfeval-template {} --threads {}".format(
		truth_vcf, query_vcf, confident_regions, regions_file, reference_fasta, output_prefix, rtg_path, rtg_template, threads)

	os.system(run_happy)

def parse_output():
	lines = open(slurm_output_file, "r").read().splitlines()

	current_sample = None
	current_platform = None
	i = 0

	while i < len(lines):
		line = lines[i]
		
		if line.startswith("Benchmarking") and "!" in line:
			current_sample = line.split()[1]
			current_platform = line.split()[2].strip("()!")
			print("Processing sample: {}, {}".format(current_sample, current_platform))
			i += 1
			continue

		if line.startswith("Benchmarking Summary:"):
			if lines[i+2].split()[0] == "INDEL" or lines[i+2].split()[0] == "SNP":
				indel_line = None
				snp_line = None

				if lines[i+2].split()[0] == "INDEL":
					indel_line = lines[i + 3]
					snp_line = lines[i + 5]
				elif lines[i + 2].split()[0] == "SNP":
					snp_line = lines[i + 3]

			if snp_line:
				snp_fields = snp_line.split()
				snp_TRUTH_TOTAL = snp_fields[2]
				snp_TRUTH_TP = snp_fields[3]
				snp_TRUTH_FN = snp_fields[4]
				snp_QUERY_FP = snp_fields[6]
				snp_RECALL = snp_fields[10]
				snp_PRECISION = snp_fields[11]
				snp_F1 = snp_fields[13]
				print("SNP F1: {}".format(snp_F1))

			if indel_line:
				indel_fields = indel_line.split()
				indel_TRUTH_TOTAL = indel_fields[2]
				indel_TRUTH_TP = indel_fields[3]
				indel_TRUTH_FN = indel_fields[4]
				indel_QUERY_FP = indel_fields[6]
				indel_RECALL = indel_fields[10]
				indel_PRECISION = indel_fields[11]
				indel_F1 = indel_fields[13]
				print("INDEL F1: {}".format(indel_F1))

			# Assign a single list instead of appending
			concordance_dict[current_sample][current_platform]["snp"] = [snp_TRUTH_TOTAL, snp_TRUTH_TP, snp_TRUTH_FN, snp_QUERY_FP, snp_RECALL, snp_PRECISION, snp_F1]
			if indel_line:
				concordance_dict[current_sample][current_platform]["indel"] = [indel_TRUTH_TOTAL, indel_TRUTH_TP, indel_TRUTH_FN, indel_QUERY_FP, indel_RECALL, indel_PRECISION, indel_F1]
			elif not indel_line:
				concordance_dict[current_sample][current_platform]["indel"] = ["NA"]*7

			if indel_line:
				i+= 6
			else:
				i+=4
			continue

		i += 1

	return concordance_dict

def write_results():
	# Define CSV headers
	header = ["sample", "platform", "type", "truth_total", "truth_tp", "truth_fn", "query_fp", "recall", "precision", "f1"]

	# Open the file for writing
	with open(parsed_results_file, mode='w', newline='') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(header)

		# Iterate over the nested dictionary
		for sample, platforms in concordance_dict.items():
			for platform, types in platforms.items():
				for type, metrics in types.items():
					row = [sample, platform, type] + metrics
					writer.writerow(row)

def main():
	for platform in platforms:
		for sample in samples:
			run_happy(platform, sample)

	parse_output()
	write_results()

if __name__ == "__main__":
	main()