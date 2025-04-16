from Bio import SeqIO
import sys

g_group_file = "hla_nom_g.txt"
hla_seq_file = "/Users/matt/Downloads/hla_nuc.fasta"

# {"A*": {'01:01:01G': [allele_1, allele_2]}}
g_group_allele_dict = dict()

# {'A*01:01:01:01': seq}
hla_seq_dict = dict()

# {'01:01:01G': seq}
g_group_seq_dict = dict()

hla_a_exon_2 = slice(73, 343)
hla_a_exon_3 = slice(343, 619)

def create_g_group_dict():
	with open(g_group_file, "r") as f:
		lines = f.read().splitlines()

	for line in lines:
		if not line.startswith("#"):
			
			fields = line.split(";")
			
			# Existing G Group Information
			if not fields[2] == "":
				gene = fields[0]
				alleles = fields[1].split("/")
				g_group = fields[2]

				if not gene in g_group_allele_dict:
					g_group_allele_dict[gene] = dict()

				g_group_allele_dict[gene][g_group] = alleles

def create_g_group_seq_dict():
	problem_counter = 0

	for record in SeqIO.parse(hla_seq_file, "fasta"):
		name = record.description.split()[1]
		seq = str(record.seq)
		hla_seq_dict[name] = seq

	for gene, g_groups in g_group_allele_dict.items():
		if gene == "B*":
			# print(f"Parsing Gene {gene}!")
			for g_group, alleles in g_groups.items():
				# print(f"Parsing G-Group: {g_group}!")
				exon_2_seq = set()
				exon_3_seq = set()
				for allele in alleles:
					allele_name = f"{gene}{allele}"
					if allele_name in hla_seq_dict:
						#if len(hla_seq_dict[allele_name]) == 1098:
						allele_seq = hla_seq_dict[allele_name]
						exon_2_seq.add(allele_seq[hla_a_exon_2])
						exon_3_seq.add(allele_seq[hla_a_exon_3])

					elif hla_seq_dict[allele_name].startswith("GCT") and hla_seq_dict[allele_name].endswith("CGG"):
						allele_seq = hla_seq_dict[allele_name]
						exon_2_seq = add(allele_seq[0:270])
						exon_3_seq = add(allele_seq[270:546])

					else:
						print(f"Allele {allele_name}, length: {len(hla_seq_dict[allele_name])}")
						
					if len(exon_2_seq) == 1 and len(exon_3_seq) == 1:
						g_group_seq_dict[g_group] = list(exon_2_seq)[0] + list(exon_3_seq)[0]
					
					else:
						print(f"Seqs not identical for all alleles in G Group {g_group}!")
						problem_counter += 1

	print(problem_counter)

def main():
	create_g_group_dict()
	create_g_group_seq_dict()
	# for gene, g_groups in g_group_allele_dict.items():
	# 	if gene == "A*":
	# 		# print(gene)
	# 		for g_group, alleles in g_groups.items():
	# 			if not g_group in g_group_seq_dict:
	# 				print(g_group)
	print(g_group_seq_dict)


if __name__ == "__main__":
	main()
