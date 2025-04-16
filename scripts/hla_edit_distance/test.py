import json
from Bio import SeqIO
import edlib
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

reference_file = "/Users/matt/Downloads/hla_nuc.fasta"
query_file = "/Users/matt/Desktop/fasta_dict.json"

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

def match_query_sequence(gene, query_name, query_sequence, alleles, n=3):
	print(f"Comparing {query_name} to {len(alleles)} HLA-{gene} alleles in {reference_file}")
	distances = []

	for allele_name, reference_sequence in alleles:
		if len(query_sequence) <= len(reference_sequence):
			dist = edlib.align(query_sequence, reference_sequence, mode="HW", task="path")["editDistance"]
		else:
			dist = edlib.align(reference_sequence, query_sequence, mode="HW", task="path")["editDistance"]
		uncertainty = dist / len(query_sequence)
		distances.append((allele_name, dist, uncertainty))

	if distances:
		top_matches = sorted(distances, key=lambda x: x[1])[:n]
		print(f"Top {n} Matches:")
		for rank, (match, dist, uncertainty) in enumerate(top_matches, start=1):
			print(f" {rank}. {match} (distance: {dist}, uncertainty: {uncertainty:.3f})")
	else:
		print(f"Error! No matches found for {query_name}")

def run_match_task(args):
	gene_abr, name, seq, allele_subset = args
	match_query_sequence(gene_abr, name, seq, allele_subset)

def process_fasta_parallel(query_file, allele_db):
	tasks = []

	with open(query_file) as f:
		fasta_data = json.load(f)

	for platform, genes in fasta_data.items():
		for gene, samples in genes.items():
			gene_abr = gene.split("-")[1]
			alleles = allele_db.get(gene_abr)

			for sample, haplotypes in samples.items():
				hap1_name = f"{sample}_HLA_{gene}_{platform}_1"
				hap2_name = f"{sample}_HLA_{gene}_{platform}_2"
				hap1_seq = haplotypes[0]
				hap2_seq = haplotypes[1]

				tasks.append((gene_abr, hap1_name, hap1_seq, alleles))
				tasks.append((gene_abr, hap2_name, hap2_seq, alleles))

	with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
		executor.map(run_match_task, tasks)

def main():
	allele_db = build_allele_database(reference_file)
	process_fasta_parallel(query_file, allele_db)

if __name__ == "__main__":
	main()
