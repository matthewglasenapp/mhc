import os
import csv
import json
import pysam
from joblib import Parallel, delayed

haploblock_dir = "/hb/scratch/mglasena/test_pacbio/processed_data/hiphase_haploblocks/"
#genes_bed = "hla_captured_genes.bed"
genes_bed = "test.bed"
vcf_dir = "/hb/scratch/mglasena/test_pacbio/processed_data/phased_vcf_hiphase/"

# "HG01891"
#samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", "IHW09409", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]
samples = ["HG002"]
genes_of_interest = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DPA1", "HLA-DPB1")

# {"gene_name": [start, stop]}
genes_dict = dict()

# {"sample_id": [[haploblock_1_start, haploblock_1_stop]]}
haploblock_dict = {sample: [] for sample in samples}

# {"sample_id": [phased_genes]}
gene_haploblock_dict = {sample: [] for sample in samples}

hla_start = 29722774
hla_stop = 33129084

# Output files
phased_genes_by_sample_csv = "phased_genes.tsv"
phased_genes_by_sample_json = "phased_genes.json"
phase_map_csv = "phase_map.csv"

# Populate dictionary of captured genes with gene name and start and stop coordinates
def create_genes_dict():
	genes = open(genes_bed, "r").read().splitlines()
	for line in genes:
		fields = line.split("\t")
		name = fields[3].split("_")[0]
		start = int(fields[1])
		stop = int(fields[2])
		genes_dict[name] = [start, stop]

# Load heterozygous variants from HiPhase VCF
def load_heterozygous_variants():
	heterozygous_sites = {sample: {"chr6": []} for sample in samples}

	for sample in samples:
		vcf_file = os.path.join(vcf_dir, f"{sample}.dedup.trimmed.hg38.chr6.phased.vcf.gz")

		if not os.path.exists(vcf_file):
			print(f"❌ Warning: VCF file missing for {sample}")
			continue  # Skip missing files

		print(f"✅ Found VCF file for {sample}")

		vcf = pysam.VariantFile(vcf_file)

		for record in vcf:
			if record.chrom != "chr6":
				continue

			genotype = record.samples[sample]["GT"]

			if genotype in [(0, 1), (1, 0)]:
				heterozygous_sites[sample]["chr6"].append(record.pos)

		print(f"➡ Loaded {len(heterozygous_sites[sample]['chr6'])} heterozygous sites for {sample}")

	print("Final heterozygous_sites dictionary keys:", heterozygous_sites.keys())  # Check which samples were added
	return heterozygous_sites

# Get list of HiPhase haploblock intervals for chromosome 6
# An alternative approach (from HiPhase): bedtools intersect -a genes_of_interest.bed -b sample.hiphase.haploblocks.bed -f 1.0 -wa 
def parse_haploblocks(sample, het_sites):
	print(f"Processing {sample}: Received {len(het_sites)} heterozygous sites")
	haploblock_list = []
	
	haploblock_file = os.path.join(haploblock_dir, f"{sample}.phased.blocks.txt")

	if not os.path.exists(haploblock_file):
		print(f"❌ Missing haploblock file for {sample}")
		return sample, []

	print(f"✅ Processing haploblock file for {sample}")

	with open(haploblock_file, "r") as f:
		haploblocks = f.read().splitlines()

	if not het_sites:
		print(f"⚠ Warning: No heterozygous sites found for {sample}, skipping haploblock extension.")

	for line in haploblocks[1:]:
		fields = line.split("\t")
		chromosome = fields[3]
		start = int(fields[4]) - 1
		stop = int(fields[5])

		if chromosome == "chr6" and stop > hla_start:
			prev_het = max([h for h in het_sites if h < start], default=start)
			next_het = min([h for h in het_sites if h > stop], default=stop)
			haploblock_list.append([prev_het - 1, next_het])

	output_bed = sample + ".haploblocks.bed"
	with open(output_bed, "w") as f:
		for haploblock in haploblock_list:
			start = haploblock[0]
			stop = haploblock[1]
			f.write(f"chr6\t{start}\t{stop}\n")
	return sample, haploblock_list


# Check whether each captured MHC gene is completely spanned by a haploblock
def evaluate_gene_haploblocks(sample, het_sites):
	gene_list = []
	haploblocks = haploblock_dict[sample]

	for gene in genes_dict:
		gene_start = genes_dict[gene][0]
		gene_stop = genes_dict[gene][1]

		gene_het_sites = [site for site in het_sites if gene_start <= site <= gene_stop]

		if not gene_het_sites:
			gene_list.append(gene)
			
			if gene in genes_of_interest:
				print(f"{gene} has no heterozygous sites in sample {sample}")
				print(f"Counting {sample} {gene} as fully phased")
			
			continue

		for block_start, block_stop in haploblocks:
			if block_start <= gene_start and block_stop >= gene_stop:
				gene_list.append(gene)
				break

	return sample, gene_list

def write_results():
	with open(phased_genes_by_sample_csv, "w", newline="") as csv_file:
		writer = csv.writer(csv_file, delimiter="\t")
		writer.writerow(["sample", "num_genes", "genes"])
		for sample, gene_list in gene_haploblock_dict.items():
			writer.writerow([sample, len(gene_list), ",".join(gene_list)])

	with open(phased_genes_by_sample_json, "w") as json_file:
		json.dump(gene_haploblock_dict, json_file, indent=4)

def make_heatmap_data():
	phased_count_dict = {gene: 0 for gene in genes_of_interest}
	phased_status_dict = {sample: [] for sample in samples}

	with open(phased_genes_by_sample_csv, "r") as f:
		records = f.read().splitlines()[1:]

	for item in records:
		fields = item.split("\t")
		sample = fields[0]
		gene_list = fields[2].split(",")

		for gene in genes_of_interest:
			if gene in gene_list:
				phased_count_dict[gene] += 1
				phased_status_dict[sample].append(1)
			else:
				phased_status_dict[sample].append(0)

	with open(phase_map_csv, "w", newline="") as csv_file:
		writer = csv.writer(csv_file, delimiter=",")
		header = ["sample", "HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DPA1", "HLA-DPB1"]
		writer.writerow(header)
		for sample, genes in phased_status_dict.items():
			writer.writerow([sample] + genes)

	prop_phased_dict = {gene: round(count / len(samples), 2) for gene, count in phased_count_dict.items()}
	for gene, prop_phased in prop_phased_dict.items():
		print(f"{gene}: {prop_phased}")



def main():
	create_genes_dict()
	heterozygous_sites = load_heterozygous_variants()

	if "HG002" not in heterozygous_sites:
		print("❌ Error: HG002 is missing from heterozygous_sites before parallel execution!")
		print("Existing keys:", heterozygous_sites.keys())
		return  # Exit early

	print("✅ HG002 exists in heterozygous_sites before parallel execution!")

	haploblocks_by_sample = Parallel(n_jobs=10)(
		delayed(parse_haploblocks)(sample, heterozygous_sites.get(sample, {}).get("chr6", [])) for sample in samples
	)

	for sample, haploblock_list in haploblocks_by_sample:
		haploblock_dict[sample] = haploblock_list


	genes_by_haploblock = Parallel(n_jobs=10)(delayed(evaluate_gene_haploblocks)(sample, heterozygous_sites.get(sample, {}).get("chr6", [])) for sample in samples)

	for sample, gene_list in genes_by_haploblock:
		gene_haploblock_dict[sample] = gene_list

	write_results()
	make_heatmap_data()
	
if __name__ == "__main__":
	main()
