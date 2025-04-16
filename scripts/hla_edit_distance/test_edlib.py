import json
from Bio import SeqIO
import edlib

reference_file = "/Users/matt/Downloads/hla_nuc.fasta"

# IHW09117 HLA-B Test Case
# query_file = "/Users/matt/Downloads/test.fa"
# Gene for Query file (e.g., "A")
gene = "B"

# All PacBio and ONT Class I haplotypes
query_file = "/Users/matt/Desktop/fasta_dict.json"

class SequenceMatcher:

	def __init__(self, reference_file):
		self.reference_file = reference_file
		self.allele_db = dict()

	def build_allele_database(self):
		for record in SeqIO.parse(self.reference_file, "fasta"):
			allele_name = record.description.split()[1]
			gene = allele_name.split("*")[0]
			sequence = str(record.seq)

			if gene not in self.allele_db:
				self.allele_db[gene] = []

			self.allele_db[gene].append((allele_name, sequence))

	def match_query_sequence(self, gene, query_name, query_sequence, n=3):
		print(f"Query Sequence Name: {query_name}")
		print(f"Query Sequence Gene: HLA-{gene}")
		alleles = self.allele_db.get(gene)
		allele_count = len(alleles)
		print(f"Comparing {query_name} to {allele_count} HLA-{gene} alleles in {reference_file}")

		distances = []

		for allele_name, reference_sequence in alleles:
			if len(query_sequence) <= len(reference_sequence):
				dist = edlib.align(query_sequence, reference_sequence, mode="HW", task="path")["editDistance"]
			elif len(query_sequence) > len(reference_sequence):
				dist = edlib.align(reference_sequence, query_sequence, mode="HW", task="path")["editDistance"]
			else:
				print("Error")
				sys.exit(1)
			uncertainty = dist / len(query_sequence)
			distances.append((allele_name, dist, uncertainty))

		top_matches = sorted(distances, key=lambda x: x[1])[:n]

		if top_matches:
			print(f"Top {n} Matches:")
			for rank, (match, dist, uncertainty) in enumerate(top_matches, start =1):
				print(f" {rank}. {match} (distance: {dist}, uncertainty: {uncertainty:.3f})")
		else:
			print(f"Error! No matches found for {query_name}")

def process_fasta(query_file, gene, matcher):
	# Input File is fasta file
	if query_file.endswith(".fasta") or query_file.endswith(".fa"):
		for record in SeqIO.parse(query_file, "fasta"):
			name = record.description
			seq = str(record.seq)
			matcher.match_query_sequence(gene, name, seq)

	# Input file is json object {platform:{gene:{samples:[alleles]}}}
	elif query_file.endswith(".json"):
		with open(query_file) as f:
			fasta_data = json.load(f)

		for platform, genes in fasta_data.items():
			for gene, samples in genes.items():
				gene_abr = gene.split("-")[1]
				for sample, haplotypes in samples.items():
					hap1_name = f"{sample}_HLA_{gene}_{platform}_1"
					hap1_seq = haplotypes[0]
					hap2_name = f"{sample}_HLA_{gene}_{platform}_2"
					hap2_seq = haplotypes[1]
					matcher.match_query_sequence(gene_abr, hap1_name, hap1_seq)
					matcher.match_query_sequence(gene_abr, hap2_name, hap2_seq)

def main():
	matcher = SequenceMatcher(reference_file)
	matcher.build_allele_database()
	process_fasta(query_file, gene, matcher)

if __name__ == "__main__":
	main()
