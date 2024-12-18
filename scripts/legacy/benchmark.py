import os

threads = 8

#samples = ["HG002", "HG003", "HG004", "HG005"]
samples = ["HG002"]

#platform = "Revio"
platform = "PromethION"

root_dir = "/hb/home/mglasena/software/happy/"
#input_vcf_dir = "/hb/scratch/ogarci12/deepvariant/old2/" + platform + "_VCF/"
#input_vcf_dir = "/hb/scratch/ogarci12/deepvariant/Revio_VCF/"
#input_vcf_dir = "/hb/scratch/mglasena/deepvariant_whathap/PromethION_VCF/HG003/"
input_vcf_dir = "/hb/scratch/mglasena/deepvariant_whahap/Revio_VCF/HG002/"

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
	#query_vcf = input_vcf_dir + sample + ".hg38.Pacbio.vcf.gz"
	#query_vcf = input_vcf_dir + sample + ".hg38.promethION.vcf.gz"
	query_vcf = input_vcf_dir + "merge_output.vcf.gz"
	truth_vcf = root_dir + sample + "/" + sample + "_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
	#truth_vcf_mhc = root_dir + "HG002/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.vcf.gz"

	if sample == "HG002" or sample == "HG003" or sample == "HG004":
		confident_regions = root_dir + sample + "/" + sample + "_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
	else:
		confident_regions = root_dir + sample + "/" + sample + "_GRCh38_1_22_v4.2.1_benchmark.bed"
	#confident_regions_MHC = root_dir + "HG002/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_callablemultinter_gt0.bed"

	run_happy = "hap.py {} {} -f {} -R {} -r {} -o {} --engine vcfeval --engine-vcfeval-path {} --engine-vcfeval-template {} --threads {}".format(
		truth_vcf, query_vcf, confident_regions, regions_file, reference_fasta, output_prefix, rtg_path, rtg_template, threads
	)

	os.system(run_happy)

def parse_output():
	pass

def main():
	for sample in samples:
		run_happy(sample)

	#parse_output()

if __name__ == "__main__":
	main()