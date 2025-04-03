import os
import json
import subprocess
import pysam

vcf2fasta = "/hb/scratch/mglasena/vcf2fasta/vcf2fasta.py"

reference_genome = "/hb/scratch/mglasena/test_pacbio/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

phased_genes = "/hb/scratch/mglasena/test_minimap/phase_results/phased_genes.whatshap.json"
phased_vcf_dir = "/hb/scratch/mglasena/test_minimap/processed_data/phased_vcf_whatshap/"
gff_dir = "/hb/scratch/mglasena/test_vcf2fasta/gff/"

samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", "IHW09409", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]

#genes_of_interest = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DPA1", "HLA-DPB1")

genes_of_interest = ["HLA-A", "HLA-B", "HLA-C"]

input_dir = os.getcwd()

filtered_vcf_dir = os.path.join(input_dir, "filtered_vcf")
os.makedirs(filtered_vcf_dir, exist_ok=True)

fasta_dir = os.path.join(input_dir, "fasta_sequences")
os.makedirs(fasta_dir, exist_ok=True)

fasta_dict = {gene: {sample: [] for sample in samples} for gene in genes_of_interest}

DNA_bases = {"A", "T", "G", "C"}

stop_codons = ["TAA", "TAG", "TGA"]

output_file = "fasta_dict.json"

def get_gff_files():
	find_cmd = "find {input_dir} -type f > gff_files.txt".format(input_dir = gff_dir)
	subprocess.run(find_cmd, shell = True, check = True)
	gff_files = open("gff_files.txt", "r").read().splitlines()
	os.remove("gff_files.txt")
	return gff_files

# Run only once
def sort_cds(gff_file):
	meta_lines = []
	cds_lines = []

	with open(gff_file, "r") as f:
		for line in f:
			if line.startswith("#"):
				meta_lines.append(line)
			else:
				fields = line.strip().split("\t")
				if fields[0] == "6":
					fields[0] = "chr6"
				if fields[2] == "CDS":
					start = int(fields[3])
					strand = fields[6]
					cds_lines.append((start, fields))

	if strand == "-":
		sorted_cds = sorted(cds_lines, key=lambda line: line[0], reverse=True)
	else: 
		sorted_cds = sorted(cds_lines, key=lambda line: line[0])

	outfile = gff_file.replace(".gff3", "_cds_sorted.gff3")

	with open(outfile, "w") as out:
		for line in meta_lines:
			out.write(line)
		for line in sorted_cds:
			fields = line[1]
			out.write("\t".join(fields) + "\n")

	print(f"Wrote: {outfile}")

def filter_vcf(sample):
	phased_vcf = os.path.join(phased_vcf_dir, sample + ".dedup.trimmed.hg38.chr6.phased.vcf.gz")
	filtered_vcf = os.path.join(filtered_vcf_dir, sample + "_filtered.vcf.gz")

	# Find the first phased variant using pysam
	vcf = pysam.VariantFile(phased_vcf)
	first_phased_pos = None
	chrom = None

	for record in vcf:
		for sample_data in record.samples.values():
			if sample_data.phased:
				first_phased_pos = record.pos
				chrom = record.chrom
				break
		if first_phased_pos:
			break

	region = f"{chrom}:{first_phased_pos}-"

	# Step 2: Filter from that region onward
	slice_cmd = "bcftools view -f PASS -r {region} {phased_vcf} -Oz -o {filtered_vcf}".format(region = region, phased_vcf = phased_vcf, filtered_vcf = filtered_vcf)
	index_cmd = "bcftools index {filtered_vcf}".format(filtered_vcf = filtered_vcf)

	subprocess.run(slice_cmd, shell=True, check=True)
	subprocess.run(index_cmd, shell=True, check=True)

def load_phased_genes():
	with open(phased_genes, "r") as f:
		data = json.load(f)

	return data

def run_vcf2fasta(sample, gene):
	gene_id = gene.lower().replace("-", "_")
	input_vcf = os.path.join(filtered_vcf_dir, sample + "_filtered.vcf.gz")
	input_gff = os.path.join(gff_dir, gene_id + "_cds_sorted.gff3")
	output_dir = os.path.join(fasta_dir, gene_id)
	os.makedirs(output_dir, exist_ok=True)
	full_output_dir = os.path.join(output_dir, sample)

	vcf2fasta_cmd = "python3 {vcf2fasta} --fasta {reference_fasta} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat CDS --blend".format(vcf2fasta = vcf2fasta, reference_fasta = reference_genome, input_vcf = input_vcf, input_gff = input_gff, output_dir = full_output_dir)
	
	subprocess.run(vcf2fasta_cmd, shell = True, check = True)

def parse_fastas():
	find_cmd = "find {input_dir} -type f > fasta_files.txt".format(input_dir = fasta_dir)
	subprocess.run(find_cmd, shell = True, check = True)
	fasta_files = open("fasta_files.txt", "r").read().splitlines()
	os.remove("fasta_files.txt")

	for file in fasta_files:
		sample = file.split("/")[-2].split("_")[0]
		gene = file.split("/")[-3].upper().replace("_", "-")
		with open(file, "r") as f:
			lines = f.read().split(">")
		allele_1 = lines[1].split("\n")[1].strip().replace("-","").strip()
		allele_2 = lines[2].split("\n")[1].strip().replace("-","").strip()

		# Remove deletion characters

		if allele_1[0:3] != "ATG" or allele_2[0:3] != "ATG":
			print("File {} does not begin with start codon!".format(file))
		
		elif not allele_1[-3:] in stop_codons or not allele_2[-3:] in stop_codons:
			print("File {} does not end with stop codon!".format(file))

		if not set(allele_1).issubset(DNA_bases):
			print("{} has invalid characters!".format(file))

		if not set(allele_2).issubset(DNA_bases):
			print("{} has invalid characters!".format(file))

		fasta_dict[gene][sample].append(allele_1)
		fasta_dict[gene][sample].append(allele_2)

def write_fasta_dict():
	with open(output_file, "w") as f:
		json.dump(fasta_dict, f, indent=2)

	print(f"Wrote fasta_dict to {output_file}")

def main():
	#for gff_file in gff_files:
		#sort_cds(gff_file)

	phased_genes = load_phased_genes()
	gff_files = get_gff_files()

	for sample in samples:
		filter_vcf(sample)

	for sample, genes in phased_genes.items():
		for gene in genes:
			if gene in genes_of_interest:
				run_vcf2fasta(sample, gene)

	parse_fastas()

	write_fasta_dict()

if __name__ == "__main__":
	main()
