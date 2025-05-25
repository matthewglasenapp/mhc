import os
import json
import subprocess
import pysam

vcf2fasta_script = "/hb/scratch/mglasena/vcf2fasta/vcf2fasta.py"
reference_genome = "/hb/groups/cornejo_lab/matt/pacbio_capture/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
gff_dir = "/hb/groups/cornejo_lab/matt/pacbio_capture/hla_gff/"
# feature = "gene"
features = ["CDS", "gene"]

samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", "IHW09409", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]

#genes_of_interest = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DPA1", "HLA-DPB1")

genes_of_interest = ["HLA-A", "HLA-B", "HLA-C"]

output_file = "fasta_dict.json"
input_dir = os.getcwd()
base_output_dir = "/hb/groups/cornejo_lab/matt/pacbio_capture/processed_data"
filtered_vcf_dir = os.path.join(base_output_dir, "vcf2fasta_vcf")
fasta_dir = os.path.join(base_output_dir, "fasta_haplotypes")
os.makedirs(filtered_vcf_dir, exist_ok=True)
os.makedirs(fasta_dir, exist_ok=True)

fasta_dict = {}
DNA_bases = {"A", "T", "G", "C"}
stop_codons = ["TAA", "TAG", "TGA"]

config = {
	"revio_hiphase": {
		"vcf_dir": "/hb/groups/cornejo_lab/matt/pacbio_capture/processed_data/phased_vcf_hiphase/",
		"vcf_suffix": ".dedup.trimmed.hg38.chr6.phased.joint.vcf.gz",
		"phased_genes": "/hb/scratch/mglasena/delete/phased_genes.hiphase.json"
	}
	# "revio_longphase": {
	# 	"vcf_dir": "/hb/groups/cornejo_lab/matt/pacbio_capture/processed_data/phased_vcf_longphase/",
	# 	"vcf_suffix": ".dedup.trimmed.hg38.chr6.phased.merged.vcf.gz",
	# 	"phased_genes": "/hb/scratch/mglasena/delete/phased_genes.longphase.json"
	# },
	# "revio_whatshap": {
	# 	"vcf_dir": "/hb/scratch/mglasena/test_minimap/processed_data/phased_vcf_whatshap/",
	# 	"vcf_suffix": ".dedup.trimmed.hg38.chr6.phased.vcf.gz",
	# 	"phased_genes": "/hb/scratch/mglasena/delete/phased_genes.whatshap.json"
	# }
	# "promethion_longphase": {
	#     "vcf_dir": "/hb/scratch/mglasena/test_ont/processed_data/phased_vcf_longphase/",
	#     "vcf_suffix": ".porechop.trimmed.hg38.rmdup.chr6.longphase.vcf.gz",
	#     "phased_genes": "/hb/scratch/mglasena/delete/phased_genes.promethion.json"
	# },
	# "promethion_whatshap": {
	#     "vcf_dir": "/hb/scratch/mglasena/test_ont/processed_data/phased_vcf_whatshap/",
	#     "vcf_suffix": ".porechop.trimmed.hg38.rmdup.chr6.phased.vcf.gz",
	#     "phased_genes": "/hb/scratch/mglasena/delete/phased_genes.promethion.json"
	# }
}

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

def filter_vcf(sample, phaser, cfg):
	input_vcf = os.path.join(cfg["vcf_dir"], sample + cfg["vcf_suffix"])
	pass_vcf = os.path.join(filtered_vcf_dir, f"{phaser}_{sample}_PASS.vcf.gz")
	filtered_vcf = os.path.join(filtered_vcf_dir, f"{phaser}_{sample}_filtered.vcf.gz")

	# Step 1: Filter for PASS variants and remove unphased hets and unsupported variant types
	filter_expr = '(GT="hom" || GT~"\\|") && (TYPE="snp" || TYPE="indel" || SVTYPE="INS" || SVTYPE="DEL") && ALT!~"^<"'
	subprocess.run(f'bcftools view -f PASS -i \'{filter_expr}\' {input_vcf} -Oz -o {pass_vcf}', shell=True, check=True)

	subprocess.run(f"bcftools index {pass_vcf}", shell=True, check=True)

	# Step 2: Find first phased variant
	vcf = pysam.VariantFile(pass_vcf)
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

	if not first_phased_pos:
		raise ValueError(f"No phased variants found in {pass_vcf}")

	region = f"{chrom}:{first_phased_pos}-"

	# Step 3: Filter VCF from the first phased variant onward
	subprocess.run(
		f"bcftools view -r {region} {pass_vcf} -Oz -o {filtered_vcf}",
		shell=True, check=True
	)
	subprocess.run(f"bcftools index {filtered_vcf}", shell=True, check=True)

def run_vcf2fasta(tag, sample, gene, feature):
	gene_id = gene.lower().replace("-", "_")
	input_vcf = os.path.join(filtered_vcf_dir, f"{tag}_{sample}_filtered.vcf.gz")
	out_dir = os.path.join(fasta_dir, gene_id, tag)
	os.makedirs(out_dir, exist_ok=True)
	full_output_dir = os.path.join(out_dir, sample)
	
	if feature == "CDS":
		input_gff = os.path.join(gff_dir, gene_id + "_cds_sorted.gff3")
		vcf2fasta_cmd = "python3 {vcf2fasta} --fasta {reference_fasta} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat CDS --blend".format(vcf2fasta = vcf2fasta_script, reference_fasta = reference_genome, input_vcf = input_vcf, input_gff = input_gff, output_dir = full_output_dir)
	elif feature == "gene":
		input_gff = os.path.join(gff_dir, gene_id + "_gene.gff3")
		vcf2fasta_cmd = "python3 {vcf2fasta} --fasta {reference_fasta} --vcf {input_vcf} --gff {input_gff} -o {output_dir} --feat gene".format(
	vcf2fasta = vcf2fasta_script, reference_fasta = reference_genome, input_vcf = input_vcf, input_gff = input_gff, output_dir = full_output_dir)
	
	subprocess.run(vcf2fasta_cmd, shell = True, check = True)

def parse_fastas():
	find_cmd = "find {input_dir} -type f > fasta_files.txt".format(input_dir = fasta_dir)
	subprocess.run(find_cmd, shell = True, check = True)
	fasta_files = open("fasta_files.txt", "r").read().splitlines()
	os.remove("fasta_files.txt")

	for file in fasta_files:
		platform = file.split("/")[-3]
		feat = file.split("/")[-2].split("_")[1]
		sample = file.split("/")[-2].split("_")[0]
		gene = file.split("/")[-4].upper().replace("_", "-")
		with open(file, "r") as f:
			lines = f.read().split(">")
		# Remove deletion characters
		allele_1 = lines[1].split("\n")[1].strip().replace("-","").strip()
		allele_2 = lines[2].split("\n")[1].strip().replace("-","").strip()

		if allele_1[0:3] != "ATG" or allele_2[0:3] != "ATG":
			print("File {} does not begin with start codon!".format(file))
		
		elif not allele_1[-3:] in stop_codons or not allele_2[-3:] in stop_codons:
			print("File {} does not end with stop codon!".format(file))

		if not set(allele_1).issubset(DNA_bases):
			print("{} has invalid characters!".format(file))

		if not set(allele_2).issubset(DNA_bases):
			print("{} has invalid characters!".format(file))

		if feat not in fasta_dict:
			fasta_dict[feat] = {}
		if platform not in fasta_dict[feat]:
			fasta_dict[feat][platform] = {}
		if gene not in fasta_dict[feat][platform]:
			fasta_dict[feat][platform][gene] = {}
		if sample not in fasta_dict[feat][platform][gene]:
			fasta_dict[feat][platform][gene][sample] = []

		fasta_dict[feat][platform][gene][sample].append(allele_1)
		fasta_dict[feat][platform][gene][sample].append(allele_2)

def write_fasta_dict():
	with open(output_file, "w") as out:
		json.dump(fasta_dict, out, indent=2)
	print(f"Wrote: {output_file}")

def main():
	# gff_files = get_gff_files()

	# for gff_file in gff_files:
	#     sort_cds(gff_file)
	
	for phaser, cfg in config.items():
		with open(cfg["phased_genes"]) as f:
			phased_genes = json.load(f)

		for sample in samples:
			filter_vcf(sample, phaser, cfg)

		for sample, genes in phased_genes.items():
			for gene in genes:
				if gene in genes_of_interest:
					for feat in features:
						run_vcf2fasta(phaser, sample, gene, feat)

	parse_fastas()
	write_fasta_dict()


if __name__ == "__main__":
	main()
