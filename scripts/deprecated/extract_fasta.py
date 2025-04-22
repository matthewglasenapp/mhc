import os
import json
import subprocess
import pysam

vcf2fasta = "/hb/scratch/mglasena/vcf2fasta/vcf2fasta.py"

reference_genome = "/hb/scratch/mglasena/test_pacbio/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

revio_phased_genes = "/hb/scratch/mglasena/test_minimap/phase_results/phased_genes.whatshap.json"
revio_phased_vcf_dir = "/hb/scratch/mglasena/test_minimap/processed_data/phased_vcf_hiphase/"
promethion_phased_genes = "/hb/scratch/mglasena/delete/phased_genes.promethion.json"
promethion_phased_vcf_dir = "/hb/scratch/mglasena/test_ont/processed_data/phased_vcf_longphase/"
gff_dir = "/hb/scratch/mglasena/test_vcf2fasta/gff/"

# platforms = ["Revio", "PromethION"]
platforms = ["Revio"]

# samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", "IHW09409", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]
samples = ["HG002"]

#genes_of_interest = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DPA1", "HLA-DPB1")

genes_of_interest = ["HLA-A", "HLA-B", "HLA-C"]

feature = "gene"
# feature = "CDS"

input_dir = os.getcwd()

filtered_vcf_dir = os.path.join(input_dir, "filtered_vcf")
os.makedirs(filtered_vcf_dir, exist_ok=True)

fasta_dir = os.path.join(input_dir, "fasta_sequences")
os.makedirs(fasta_dir, exist_ok=True)

fasta_dict = {platform: {gene: {sample: [] for sample in samples} for gene in genes_of_interest} for platform in platforms}

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
	gene_line = []
	strand = None

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
				elif fields[2] == "gene":
					gene_line.append(fields)


	if strand == "-":
		sorted_cds = sorted(cds_lines, key=lambda line: line[0], reverse=True)
	else: 
		sorted_cds = sorted(cds_lines, key=lambda line: line[0])

	outfile_cds = gff_file.replace(".gff3", "_cds_sorted.gff3")
	outfile_gene = gff_file.replace(".gff3", "_gene.gff3")

	with open(outfile_cds, "w") as out:
		for line in meta_lines:
			out.write(line)
		for line in sorted_cds:
			fields = line[1]
			out.write("\t".join(fields) + "\n")

	with open(outfile_gene, "w") as out2:
		for line in meta_lines:
			out2.write(line)
		out2.write("\t".join(gene_line[0]) + "\n")

	print(f"Wrote: {outfile_gene}")
	print(f"Wrote: {outfile_cds}")

def filter_vcf(sample, platform):
	if platform == "Revio":
		#phased_vcf = os.path.join(revio_phased_vcf_dir, sample + ".dedup.trimmed.hg38.chr6.phased.vcf.gz")
		phased_vcf = os.path.join(revio_phased_vcf_dir, sample + ".dedup.trimmed.hg38.chr6.phased.merged.vcf.gz")
	elif platform == "PromethION":
		phased_vcf = os.path.join(promethion_phased_vcf_dir, sample + ".porechop.trimmed.hg38.rmdup.chr6.longphase.vcf.gz")

	pass_only_vcf = os.path.join(filtered_vcf_dir, f"{platform}_{sample}_PASS.vcf.gz")
	filtered_vcf = os.path.join(filtered_vcf_dir, f"{platform}_{sample}_filtered.vcf.gz")

	# --include 'FMT/GQ>=20' -f PASS 
	get_pass_vcf_cmd = "bcftools view -f PASS {input_vcf} -Oz -o {output_vcf}".format(input_vcf = phased_vcf, output_vcf = pass_only_vcf)
	
	index_cmd = "bcftools index {input_vcf}".format(input_vcf = pass_only_vcf)

	subprocess.run(get_pass_vcf_cmd, shell=True, check=True)
	subprocess.run(index_cmd , shell=True, check=True)

	# Find the first phased variant using pysam
	vcf = pysam.VariantFile(pass_only_vcf)
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
	slice_cmd = "bcftools view -r {region} {input_vcf} -Oz -o {output_vcf}".format(region = region, input_vcf = pass_only_vcf, output_vcf = filtered_vcf)
	index_cmd = "bcftools index {input_vcf}".format(input_vcf = filtered_vcf)

	subprocess.run(slice_cmd, shell=True, check=True)
	subprocess.run(index_cmd, shell=True, check=True)

	os.remove(pass_only_vcf)
	os.remove(pass_only_vcf + ".csi")

def load_phased_genes(platform):
	if platform == "Revio":
		phased_genes = revio_phased_genes
	elif platform == "PromethION":
		phased_genes = promethion_phased_genes

	with open(phased_genes, "r") as f:
		data = json.load(f)

	return data

def run_vcf2fasta(platform, sample, gene):
	gene_id = gene.lower().replace("-", "_")
	input_vcf = os.path.join(filtered_vcf_dir, f"{platform}_{sample}_filtered.vcf.gz")
	output_dir = os.path.join(fasta_dir, gene_id, platform)
	os.makedirs(output_dir, exist_ok=True)
	full_output_dir = os.path.join(output_dir, sample)
	
	if feature == "CDS":
		input_gff = os.path.join(gff_dir, gene_id + "_cds_sorted.gff3")
		vcf2fasta_cmd = "python3 {vcf2fasta} --fasta {reference_fasta} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat CDS --blend".format(vcf2fasta = vcf2fasta, reference_fasta = reference_genome, input_vcf = input_vcf, input_gff = input_gff, output_dir = full_output_dir)
	elif feature == "gene":
		input_gff = os.path.join(gff_dir, gene_id + "_gene.gff3")
		vcf2fasta_cmd = "python3 {vcf2fasta} --fasta {reference_fasta} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat gene".format(vcf2fasta = vcf2fasta, reference_fasta = reference_genome, input_vcf = input_vcf, input_gff = input_gff, output_dir = full_output_dir)
	
	subprocess.run(vcf2fasta_cmd, shell = True, check = True)

def parse_fastas():
	find_cmd = "find {input_dir} -type f > fasta_files.txt".format(input_dir = fasta_dir)
	subprocess.run(find_cmd, shell = True, check = True)
	fasta_files = open("fasta_files.txt", "r").read().splitlines()
	os.remove("fasta_files.txt")

	for file in fasta_files:
		platform = file.split("/")[-3]
		sample = file.split("/")[-2].split("_")[0]
		gene = file.split("/")[-4].upper().replace("_", "-")
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

		if platform not in fasta_dict:
			fasta_dict[platform] = {}
		if gene not in fasta_dict[platform]:
			fasta_dict[platform][gene] = {}
		if sample not in fasta_dict[platform][gene]:
			fasta_dict[platform][gene][sample] = []


		fasta_dict[platform][gene][sample].append(allele_1)
		fasta_dict[platform][gene][sample].append(allele_2)

def write_fasta_dict():
	with open(output_file, "w") as f:
		json.dump(fasta_dict, f, indent=2)

	print(f"Wrote fasta_dict to {output_file}")

def main():
	# gff_files = get_gff_files()

	# for gff_file in gff_files:
	# 	sort_cds(gff_file)

	for platform in platforms:
		phased_genes = load_phased_genes(platform)

		for sample in samples:
			filter_vcf(sample, platform)

		for sample, genes in phased_genes.items():
			for gene in genes:
				if gene in genes_of_interest:
					run_vcf2fasta(platform, sample, gene)

	parse_fastas()

	write_fasta_dict()

if __name__ == "__main__":
	main()
