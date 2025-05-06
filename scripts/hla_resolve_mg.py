from Bio import SeqIO
import numpy as np
from numba import njit

reference_file = "/Users/matt/Downloads/hla_gen.fasta"
input_file = "/Users/matt/Desktop/drb1.fasta"

BASE_MAP = np.full(256, 4, dtype=np.uint8)
BASE_MAP[ord('A')] = 0
BASE_MAP[ord('C')] = 1
BASE_MAP[ord('G')] = 2
BASE_MAP[ord('T')] = 3
BASE_MAP[ord('a')] = 0
BASE_MAP[ord('c')] = 1
BASE_MAP[ord('g')] = 2
BASE_MAP[ord('t')] = 3

REVERSE_MAP = np.array(['A', 'C', 'G', 'T', 'N'])

def encode_sequence(fasta_sequence):
	return BASE_MAP[np.frombuffer(fasta_sequence.encode('ascii'), dtype=np.uint8)]

def decode_sequence(encoded_array):
	return ''.join(REVERSE_MAP[encoded_array])

def check_non_acgt_bases(allele_db):
	non_acgt_count = 0
	for gene, allele_list in allele_db.items():
		for allele_name, encoded_seq in allele_list:
			if np.any(encoded_seq == 4):
				non_acgt_count += 1
				print(f"[WARNING] Non-ACGT bases found in allele: {allele_name} (gene {gene})")
	if non_acgt_count == 0:
		print("All sequences are clean — no non-ACGT bases found.")

@njit
def calculate_edit_distance(query_seq, reference_seq, max_dist):
	m, n = len(query_seq), len(reference_seq)

	if abs(m - n) > max_dist:
		return max_dist + 1

	current = np.arange(n + 1, dtype=np.int32)
	previous = np.zeros(n + 1, dtype=np.int32)

	for i in range(1, m + 1):
		previous, current = current, previous
		current[0] = i
		min_val = i

		for j in range(1, n + 1):
			if query_seq[i-1] == reference_seq[j-1]:
				current[j] = previous[j-1]
			else:
				current[j] = 1 + min(previous[j], current[j-1], previous[j-1])
			min_val = min(min_val, current[j])

		if min_val > max_dist:
			return max_dist + 1

	return current[n] if current[n] <= max_dist else max_dist + 1

class SequenceMatcher:

	def __init__(self, reference_file):
		self.reference_file = reference_file
		self.allele_db = dict()

	def build_allele_database(self):
		for record in SeqIO.parse(self.reference_file, "fasta"):
			allele_name = record.description.split()[1]
			gene = allele_name.split("*")[0]
			encoded_seq = encode_sequence(str(record.seq))

			if gene not in self.allele_db:
				self.allele_db[gene] = []

			self.allele_db[gene].append((allele_name, encoded_seq))

	def match_query_sequence(self, gene, query_name, query_sequence, max_dist=1000):
		print(f"Query Sequence Name: {query_name}")
		print(f"Query Sequence Gene: HLA-{gene}")
		alleles = self.allele_db.get(gene)
		allele_count = len(alleles)
		print(f"Comparing {query_name} to {allele_count} HLA-{gene} alleles in {input_file}")

		encoded_query = encode_sequence(query_sequence)
		best_match = None
		best_distance = float('inf')

		for allele_name, encoded_ref in alleles:
			dist = calculate_edit_distance(encoded_query, encoded_ref, max_dist)
			if dist < best_distance:
				best_distance = dist
				best_match = allele_name

		if best_match:
			print(f"Best Match: {best_match} (distance = {best_distance})")
		else:
			print("No match found within max_dist")
	
def process_fasta(input_file, gene, matcher):
	input_sequences = []

	for record in SeqIO.parse(input_file, "fasta"):
		name = record.description
		seq = str(record.seq)
		input_sequences.append((name,seq))

	for query_name, query_seq in input_sequences:
		matches = matcher.match_query_sequence(gene, query_name, query_seq)

def main():
	matcher = SequenceMatcher(reference_file)
	matcher.build_allele_database()
	check_non_acgt_bases(matcher.allele_db)
	process_fasta(input_file, "DRB1", matcher)

if __name__ == "__main__":
	main()

