import os
import csv
import json
import pysam
from joblib import Parallel, delayed

# For HiPhase
haploblock_dir = "/hb/scratch/mglasena/test_pacbio/processed_data/phased_vcf_hiphase/"
vcf_dir = "/hb/scratch/mglasena/test_pacbio/processed_data/phased_vcf_hiphase/"

# For WhatsHap
#haploblock_dir = "/hb/scratch/mglasena/test_pacbio/processed_data/phased_vcf_whatshap"
#vcf_dir = "/hb/scratch/mglasena/test_pacbio/processed_data/phased_vcf_whatshap"

# For PromethION WhatsHap
#haploblock_dir = "/hb/scratch/mglasena/MHC/scripts/haploblocks/"
#vcf_dir = "/hb/scratch/mglasena/MHC/genotypes/promethion/"

#genes_bed = "hla_captured_genes.bed"
genes_bed = "test.bed"

# "HG01891"
samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", "IHW09409", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]
genes_of_interest = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DPA1", "HLA-DPB1")

# {"gene_name": [start, stop]}
genes_dict = dict()

# {"sample_id": [[haploblock_1_start, haploblock_1_stop]]}
haploblock_dict = {sample: [] for sample in samples}

# {"sample_id": [phased_genes]}
gene_haploblock_dict = {sample: [] for sample in samples}

incomplete_data = []

# Extended MHC coordinates
mhc_start = 29555628
mhc_stop = 33409896

# Output files
phased_genes_by_sample_csv = "phased_genes.tsv"
phased_genes_by_sample_json = "phased_genes.json"
phase_map_csv = "phase_map.csv"

# For unphased genes of interest, [gene, sample, haploblocks, prop_bases_phased]
incomplete_file = "incomplete.csv"

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
		# For HiPhase/WhatsHap Revio
		vcf_file = os.path.join(vcf_dir, f"{sample}.dedup.trimmed.hg38.chr6.SNV.phased.vcf.gz")
		
		# For old PromethION WhatsHap
		#vcf_file = os.path.join(vcf_dir, sample, f"{sample}.hg38.promethion.phased.vcf.gz")

		vcf = pysam.VariantFile(vcf_file)

		for record in vcf:
			if record.chrom != "chr6":
				continue

			if record.pos < mhc_start or record.pos > mhc_stop:
				continue

			# For Revio HiPhase/WhatsHap
			genotype = record.samples[sample]["GT"]

			# For old WhatsHap
			#sample_name = list(vcf.header.samples)[0]
			#genotype = record.samples[sample_name]["GT"]

			if genotype in [(0, 1), (1, 0)]:
				heterozygous_sites[sample]["chr6"].append(record.pos)

		print(f"Sample {sample} has {len(heterozygous_sites[sample]['chr6'])} heterozygous extended MHC genotypes")

	return heterozygous_sites

# Get list of HiPhase haploblock intervals for MHC
def parse_haploblocks(sample, het_sites):
	haploblock_list = []
	
	# For HiPhase
	haploblock_file = os.path.join(haploblock_dir, f"{sample}.phased.blocks.txt")
	# For WhatsHap[]
	#haploblock_file = os.path.join(haploblock_dir, f"{sample}.phased.haploblocks.txt")
	# For PromethION WhatsHap
	#haploblock_file = os.path.join(haploblock_dir, f"{sample}_promethion_haploblocks.tsv")

	print(f"Parsing {sample} haploblock file!")

	with open(haploblock_file, "r") as f:
		haploblocks = f.read().splitlines()

	for line in haploblocks[1:]:
		fields = line.split("\t")
		#For HiPhase
		chromosome = fields[3]
		start = int(fields[4]) - 1
		stop = int(fields[5])
		# For WhatsHap
		# chromosome = fields[1]
		# start = int(fields[3]) - 1
		# stop = int(fields[4])

		if chromosome == "chr6" and stop > mhc_start:
			haploblock_list.append([start,stop])

	return sample, haploblock_list

# Check whether each captured MHC gene is completely spanned by a haploblock
def evaluate_gene_haploblocks(sample, het_sites):
	# List of fully phased genes
	gene_list = []
	
	# List of genes with partially overlapping haploblock
	sample_incomplete_data = []
	

	haploblocks = haploblock_dict[sample]

	for gene in genes_dict:
		gene_start = genes_dict[gene][0]
		gene_stop = genes_dict[gene][1]
		gene_length = gene_stop - gene_start

		gene_het_sites = [site for site in het_sites if site >= gene_start and site <= gene_stop]

		# If the gene has 0 or 1 heterozygous sites, it is effectively fully phased
		if len(gene_het_sites) <= 1:
			gene_list.append(gene)
			continue

		# If the gene is completely spanned by a single haploblock, it is fully phased 
		fully_phased = False
		for block_start, block_stop in haploblocks:
			if block_start <= gene_start and block_stop >= gene_stop:
				gene_list.append(gene)
				fully_phased = True
				break

			# Check to see if unphased genes become fully phased when extending haploblocks through homozygous regions
			# Find first heterozygous site upstream of block start. Do not extend if no heterozygous sites. 
			if block_start <= gene_stop and block_stop >= gene_start:
				extended_start = max([h for h in het_sites if h < block_start], default=block_start)
				extended_start = max(extended_start, gene_start)

				# Find first heterozgyous site downstream of block stop. Do not extend if no heterozygous sites
				extended_stop = min([h for h in het_sites if h > block_stop], default=block_stop)
				extended_stop = min(extended_stop, gene_stop)

				# Consider gene fully phased if haploblock extension through homozygous bases fully spans gene coordinates. 
				if extended_start <= gene_start and extended_stop >= gene_stop:
					gene_list.append(gene)
					fully_phased = True
					break

		# Get details on haploblock overlap for genes of interest (HLA Class I/II) that were not fully phased 
		if not fully_phased and gene in genes_of_interest:
			overlapping_haploblocks = [(block_start, block_stop) for (block_start, block_stop) in haploblocks if (block_stop >= gene_start and block_start <= gene_stop)]
			num_pre_merge_haploblocks = len(overlapping_haploblocks)  # Track count before extension & merging
			print(f"Processing {sample} {gene}")
			print(f"Overlapping unextended haploblocks: {len(overlapping_haploblocks)}")

			# Step 1: Extend each haploblock independently
			extended_haploblocks = []
			for block_start, block_stop in overlapping_haploblocks:
				extended_start = max([h for h in het_sites if h < block_start], default=block_start)
				extended_start = max(extended_start, gene_start)

				extended_stop = min([h for h in het_sites if h > block_stop], default=block_stop)
				extended_stop = min(extended_stop, gene_stop)

				extended_haploblocks.append((extended_start, extended_stop))

			print(f"Extended haploblocks: {len(extended_haploblocks)}")

			# Step 2: Merge overlapping extended haploblocks
			merged_intervals = []
			extended_haploblocks.sort()

			for start, stop in extended_haploblocks:
				if not merged_intervals or start > merged_intervals[-1][1]:
					merged_intervals.append((start, stop))
				else:
					last_start, last_stop = merged_intervals[-1]
					new_stop = max(last_stop, stop)
					merged_intervals[-1] = (last_start, new_stop)

			print(f"Merged haploblocks: {len(merged_intervals)}")

			# Step 3: Compute overlap & percentage
			total_overlap = 0
			for start, stop in merged_intervals:
				overlap_start = max(start, gene_start)
				overlap_stop = min(stop, gene_stop)
				overlap_length = max(0, overlap_stop - overlap_start)
				total_overlap += overlap_length
			prop_overlap = total_overlap / gene_length
			prop_phased_string = f"{prop_overlap * 100:.2f}%"

			# Use the pre-merge haploblock count, not merged count!
			sample_incomplete_data.append([sample, gene, num_pre_merge_haploblocks, prop_phased_string])
			print(f"{sample} {gene}, Pre-Merge Haploblocks: {num_pre_merge_haploblocks}, Prop Phased: {prop_phased_string}")

	return sample, gene_list, sample_incomplete_data

def write_results():
	# Write fully phased genes
	with open(phased_genes_by_sample_csv, "w", newline="") as csv_file:
		writer = csv.writer(csv_file, delimiter="\t")
		writer.writerow(["sample", "num_genes", "genes"])
		for sample, gene_list in gene_haploblock_dict.items():
			writer.writerow([sample, len(gene_list), ",".join(gene_list)])

	with open(phased_genes_by_sample_json, "w") as json_file:
		json.dump(gene_haploblock_dict, json_file, indent=4)

	# Write incomplete.csv if there are entries 
	if incomplete_data:
		with open(incomplete_file, "w", newline="") as csvfile:
			csv_writer = csv.writer(csvfile)
			csv_writer.writerow(["sample", "gene", "num_haploblocks", "prop_phased"])
			csv_writer.writerows(incomplete_data)

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

	haploblocks_by_sample = Parallel(n_jobs=10)(
		delayed(parse_haploblocks)(sample, heterozygous_sites.get(sample, {}).get("chr6", [])) for sample in samples)

	for sample, haploblock_list in haploblocks_by_sample:
		haploblock_dict[sample] = haploblock_list

	genes_by_haploblock = Parallel(n_jobs=10)(
		delayed(evaluate_gene_haploblocks)(sample, heterozygous_sites.get(sample, {}).get("chr6", [])) for sample in samples)

	for sample, gene_list, sample_incomplete_data in genes_by_haploblock:
		gene_haploblock_dict[sample] = gene_list
		incomplete_data.extend(sample_incomplete_data)

	write_results()
	make_heatmap_data()
	
if __name__ == "__main__":
	main()
