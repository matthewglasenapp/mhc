import csv
from Bio import SeqIO
import edlib
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

reference_file = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/IPD_IMGT/hla_nuc.fasta"
query_file = "HLA_Class_I_haplotypes.fa"
output_csv = "hla_typing_results.csv"

# Create dictionary of reference alleles in the format of {gene: [(allele_name, allele_seq)]
def build_allele_database(reference_file):
	allele_db = dict()
	for record in SeqIO.parse(reference_file, "fasta"):
		allele_name = record.description.split()[1]
		gene = allele_name.split("*")[0]
		sequence = str(record.seq)
		if gene not in allele_db:
			allele_db[gene] = []
		allele_db[gene].append((allele_name, sequence))
	return allele_db

# Match query sequence to all reference alleles for a given gene using edlib.align()
# Return the name of the best match 
def match_query_sequence(gene, query_sequence, alleles, n=3):
	distances = []

	for allele_name, reference_sequence in alleles:
		if len(query_sequence) <= len(reference_sequence):
			dist = edlib.align(query_sequence, reference_sequence, mode="HW", task="path")["editDistance"]
		else:
			dist = edlib.align(reference_sequence, query_sequence, mode="HW", task="path")["editDistance"]

		uncertainty = dist / len(query_sequence)
		distances.append((allele_name, dist, uncertainty))

	top_matches = sorted(distances, key=lambda x: x[1])[:n]
	best_match = top_matches[0][0]

	return best_match

# Worker function to match a query haplotype sequence to the best-matching reference allele.
def run_match_task(args):
	sample, tag, gene, haplotype, seq, reference_alleles = args
	best_match = match_query_sequence(gene, seq, reference_alleles)
	return (sample, tag, gene, haplotype, best_match)

# Parses a FASTA file containing haplotype sequences, comparing each sequence to all reference alleles using multiprocessing
def process_fasta_to_csv(query_file, allele_db):
	tasks = []

	for record in SeqIO.parse(query_file, "fasta"):
		query_name = record.id
		query_name_fields = record.id.split("_")
		sample = query_name_fields[0]
		gene = query_name_fields[1]
		gene_abr = gene.split("-")[1]
		platform = query_name_fields[2]
		phaser = query_name_fields[3]
		haplotype = query_name_fields[4]
		tag = f"{platform}_{phaser}"
		seq = str(record.seq)
		reference_alleles = allele_db[gene_abr]
		tasks.append((sample, tag, gene_abr, haplotype, seq, reference_alleles))

	with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
		results = list(executor.map(run_match_task, tasks))

	return results

# Organize the matching results into a dictionary for parsing
# data = {platform: {gene: {sample: [hap1_best_match, hap2_best_match]}}}
def organize_results(results):
	data = {}

	for sample, platform, gene, haplotype, best_match in results:
		if platform not in data:
			data[platform] = {}
		if gene not in data[platform]:
			data[platform][gene] = {}
		if sample not in data[platform][gene]:
			data[platform][gene][sample] = [None, None]

		hap_idx = int(haplotype) - 1
		data[platform][gene][sample][hap_idx] = best_match

	return data

# Write the best matches to CSV
def write_results_to_csv(data, output_file):
	genes = sorted({gene for platform in data for gene in data[platform]})
	samples = sorted({sample for platform in data for gene in data[platform] for sample in data[platform][gene]})

	output_columns = [f"{gene}_{i}" for gene in genes for i in (1, 2)]

	with open(output_file, "w", newline="") as f:
		writer = csv.writer(f)
		writer.writerow(["Sample", "Platform"] + output_columns)

		for platform in sorted(data):
			for sample in samples:
				row = [sample, platform]

				for gene in genes:
					haplotypes = data[platform][gene][sample]
					row.extend(haplotypes)

				writer.writerow(row)


def main():
	allele_db = build_allele_database(reference_file)
	matching_results = process_fasta_to_csv(query_file, allele_db)
	organized_data = organize_results(matching_results)
	write_results_to_csv(organized_data, output_csv)

if __name__ == "__main__":
	main()
