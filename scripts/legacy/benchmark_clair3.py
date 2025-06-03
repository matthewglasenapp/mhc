import os

threads = 8

samples = ["HG002", "HG003", "HG004", "HG005"]

platform = "PromethION"

root_dir = "/hb/home/mglasena/software/happy/"
input_vcf_dir = "/hb/scratch/mglasena/mhc_genotype_calls/PromethION_VCF/"

output_dir = "/hb/scratch/mglasena/hap_py_results"
reference_fasta = root_dir + "Homo_sapiens.GRCh38.dna.primary_assembly_renamed.fa"
regions_file = root_dir + "merged_hla_legacy.bed"

# Path to rtg tools
rtg_path = "/hb/home/mglasena/.conda/envs/happy/bin/rtg"
rtg_template = root_dir + "rtg_sdf_template"

# Create output directory if it doesn't exist
make_output_dir = "mkdir -p {}".format(output_dir)
os.system(make_output_dir)

def run_happy(sample):
	print("Benchmarking {}!".format(sample))

	output_prefix = output_dir + "/" + platform + "_" + sample
	query_vcf = input_vcf_dir + sample + "/merge_output.vcf.gz"
	truth_vcf = root_dir + sample + "/" + sample + "_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

	if sample == "HG002" or sample == "HG003" or sample == "HG004":
		confident_regions = root_dir + sample + "/" + sample + "_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
	else:
		confident_regions = root_dir + sample + "/" + sample + "_GRCh38_1_22_v4.2.1_benchmark.bed"

	run_happy = "hap.py {} {} -f {} -R {} -r {} -o {} --engine vcfeval --engine-vcfeval-path {} --engine-vcfeval-template {} --threads {}".format(
		truth_vcf, query_vcf, confident_regions, regions_file, reference_fasta, output_prefix, rtg_path, rtg_template, threads
	)

	os.system(run_happy)

def main():
	for sample in samples:
		run_happy(sample)

if __name__ == "__main__":
	main()