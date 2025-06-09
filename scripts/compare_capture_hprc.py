import edlib
import json
import sys
import csv
import os
import subprocess
from Bio import SeqIO

samples = ["HG002", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "NA19240", "NA20129", "NA21309"]
genes = ["HLA-A", "HLA-B", "HLA-C"]
platforms = ["revio_hiphase", "revio_longphase", "promethion_longphase"]
fasta_dict = "fasta_dict.json"
hprc_fasta_dir = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/hprc/"
feature = "gene"

HPRC_dict = dict()

def create_HPRC_dict():
	find_cmd = f"find {hprc_fasta_dir} -name '*.fa' -type f"
	fasta_files = subprocess.run(find_cmd, capture_output=True, text=True, shell=True).stdout.strip().split("\n")

	for file in fasta_files:
		sample = os.path.basename(file).split(".")[0]

		for record in SeqIO.parse(file, "fasta"):
			gene = record.id.strip().split("::")[0]
			seq = str(record.seq)

			if gene not in HPRC_dict:
				HPRC_dict[gene] = {}

			if sample not in HPRC_dict[gene]:
				HPRC_dict[gene][sample] = ()

			HPRC_dict[gene][sample] += (seq,)	

def create_capture_seq_dict():
	with open(fasta_dict, "r") as f:
		data = json.load(f)
	return data

def compute_edit_distance(hap1, hap2):
	if len(hap1) > len(hap2):
		results = edlib.align(hap2, hap1, mode="HW", task="path")
	else:
		results = edlib.align(hap1, hap2, mode="HW", task="path")
	
	return results

def get_edlib_alignment_stats(edlib_alignment):
	edit_distance = edlib_alignment["editDistance"]
	start, end = edlib_alignment["locations"][0]
	alignment_length = end - start + 1
	match_identity = 1 - (edit_distance / alignment_length)
	return edit_distance, alignment_length, match_identity

def compare_haplotypes(sample, gene, capture_alleles):
	if gene not in HPRC_dict or sample not in HPRC_dict[gene]:
		print(f"Missing HPRC alleles for {sample} {gene}")
		return None
	
	HPRC_alleles = HPRC_dict[gene][sample]

	if len(HPRC_alleles) == 1:
		alignment_1 = compute_edit_distance(capture_alleles[0], HPRC_alleles[0])
		alignment_2 = compute_edit_distance(capture_alleles[1], HPRC_alleles[0])
		if alignment_1["editDistance"] < alignment_2["editDistance"]:
			best_alignment = alignment_1
		elif alignment_2["editDistance"] < alignment_1["editDistance"]:
			best_alignment = alignment_2
		else:
			print("Edit Distance Identical")
			sys.exit(1)

		dist, aln_len, identity = get_edlib_alignment_stats(best_alignment)
		print((dist, aln_len, identity))
		return (dist, aln_len, identity)

	one_one = compute_edit_distance(capture_alleles[0], HPRC_alleles[0])
	two_two = compute_edit_distance(capture_alleles[1], HPRC_alleles[1])
	option_one_sum = one_one["editDistance"] + two_two["editDistance"]

	one_two = compute_edit_distance(capture_alleles[0], HPRC_alleles[1])
	two_one = compute_edit_distance(capture_alleles[1], HPRC_alleles[0])
	option_two_sum = one_two["editDistance"] + two_one["editDistance"]

	if option_one_sum < option_two_sum:
		best_pair = (one_one, two_two)
	elif option_two_sum < option_one_sum:
		best_pair = (one_two, two_one)
	else:
		print("Option 1 and Option 2 are Identical!")
		sys.exit(1)

	alignment_1_metrics = get_edlib_alignment_stats(best_pair[0])
	alignment_2_metrics = get_edlib_alignment_stats(best_pair[1])
	print(alignment_1_metrics)
	print(alignment_2_metrics)
	return (alignment_1_metrics[0], alignment_1_metrics[1], alignment_1_metrics[2], alignment_2_metrics[0], alignment_2_metrics[1], alignment_2_metrics[2])


def main():
	create_HPRC_dict()
	capture_seq_dict = create_capture_seq_dict()
	
	with open("allele_match.csv", "w") as f:
		writer = csv.writer(f)
		writer.writerow(["platform", "gene", "sample", "hap_1_dist", "hap_2_dist", "alignment_length_1", "alignment_length_2", "hap_1_id", "hap_2_id"])

		for platform in platforms:
			for gene in genes:
				for sample in samples:
					capture_data = capture_seq_dict[feature][platform][gene][sample]
					result = compare_haplotypes(sample, gene, capture_data)

					if result is None:
						continue

					if len(HPRC_dict[gene][sample]) == 1:
						hap_1_dist, hap_1_len, hap_1_id = result
						hap_2_dist, hap_2_len, hap_2_id = "NA", "NA", "NA"
					else:
						hap_1_dist, hap_1_len, hap_1_id, hap_2_dist, hap_2_len, hap_2_id = result

					writer.writerow([platform, gene, sample, hap_1_dist, hap_2_dist, hap_1_len, hap_2_len, hap_1_id, hap_2_id])

if __name__ == "__main__":
	main()
