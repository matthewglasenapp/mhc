ihw_data_file = "/Users/matt/Downloads/IHW_samples.txt"
ihw_seq_file = "/Users/matt/Downloads/hla_nuc.fasta"

ihw_data_dict = dict()
ihw_seq_dict = dict()

allele_database_counts = {
	"A": 0,
	"B": 0,
	"C": 0,
	"DPA1": 0,
	"DPB1": 0,
	"DQA1": 0,
	"DQB1": 0,
	"DRB1": 0,
	"DRB3": 0,
	"DRB4": 0,
	"DRB5": 0
}

def populate_ihw_dict():
	with open(ihw_data_file, "r") as f:
		lines = f.read().splitlines()

	for line in lines:
		fields = line.split("\t")
		fields = [f if f.strip() else "NA" for f in fields]
		availability, ihw_id, alternate_id, sex, ethnicity, workshop = fields[:6]
		A_1, A_2, B_1, B_2, C_1, C_2, DPA1_1, DPA1_2, DPB1_1, DPB1_2, DQA1_1, DQA1_2, DQB1_1, DQB1_2, DRB1_1, DRB1_2, DRB3_1, DRB3_2, DRB4_1, DRB4_2, DRB5_1, DRB5_2 = fields[6:29]
		
		ihw_data_dict[ihw_id] = {
			"availability": availability,
			"alternate_id": alternate_id,
			"sex": sex,
			"ethnicity": ethnicity,
			"workshop": workshop,
			"hla": {
				"A": [A_1, A_2],
				"B": [B_1, B_2],
				"C": [C_1, C_2],
				"DPA1": [DPA1_1, DPA1_2],
				"DPB1": [DPB1_1, DPB1_2],
				"DQA1": [DQA1_1, DQA1_2],
				"DQB1": [DQB1_1, DQB1_2],
				"DRB1": [DRB1_1, DRB1_2],
				"DRB3": [DRB3_1, DRB3_2],
				"DRB4": [DRB4_1, DRB4_2],
				"DRB5": [DRB5_1, DRB5_2]
			}
		}

def populate_seq_dict():
	with open(ihw_seq_file, "r") as f:
		seqs = f.read().split(">")
	for seq in seqs[1:]:
		info = seq.split("\n",1)[0].split(" ")
		ID = info[0]
		name = info[1]
		length = info[2]
		seq = seq.split("\n",1)[1].replace("\n", "").strip()
		ihw_seq_dict[name] = seq

def count_seqs():
	alleles = list(ihw_seq_dict.keys())
	for allele in alleles:
		gene = allele.split("*")[0]
		if gene in allele_database_counts:
			allele_database_counts[gene] += 1
		# else:
		# 	print(allele)
		# 	print(gene)

def main():
	populate_ihw_dict()
	populate_seq_dict()
	count_seqs()

	# print(len(ihw_seq_dict))
	print(allele_database_counts)

	IHW09117_allele_1_id = ihw_data_dict["IHW09117"]["hla"]["C"][0]
	IHW09117_allele_2_id = ihw_data_dict["IHW09117"]["hla"]["C"][1]

	IHW09117_allele_1_seq = ihw_seq_dict[ihw_data_dict["IHW09117"]["hla"]["C"][0]]
	IHW09117_allele_2_seq = ihw_seq_dict[ihw_data_dict["IHW09117"]["hla"]["C"][1]]

	print(IHW09117_allele_1_id)
	print(IHW09117_allele_2_id)
	# print(IHW09117_allele_1_seq)
	# print(IHW09117_allele_2_seq)



if __name__ == "__main__":
	main()

