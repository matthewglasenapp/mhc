import csv

spechla_pacbio = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/spechla_results/PacBio_spechla_results.txt"
spechla_ont = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/spechla_results/ONT_spechla_results.txt"
hla_resolve = "hla_typing_results.csv"
ihw_data_file = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/IPD_IMGT/IHW_samples.txt"
output_file = "typing_comparison.csv"

spechla_platforms = ["revio", "promethion"]
resolve_platforms = ["revio_hiphase", "revio_longphase", "promethion_longphase"]

genes = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", "IHW09409", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]

# Initialize with separate dictionaries
spechla = {p: {g: {s: [] for s in samples} for g in genes} for p in spechla_platforms}
resolve = {p: {g: {s: [] for s in samples} for g in genes} for p in resolve_platforms}

ihw_data_dict = dict()

def populate_spechla_dict():
	with open(spechla_pacbio, "r") as f:
		pacbio_lines = f.read().splitlines()

	for line in pacbio_lines[2:]:
		platform = "revio"
		sample = line.split("\t")[0]
		data = line.split("\t")[1:]
		for item in data:
			gene = item.split("*")[0]
			allele = item
			if gene in spechla[platform] and sample in spechla[platform][gene]:
				spechla[platform][gene][sample].append(allele)

	with open(spechla_ont, "r") as f:
		ont_lines = f.read().splitlines()

	for line in ont_lines[2:]:
		platform = "promethion"
		sample = line.split("\t")[0]
		data = line.split("\t")[1:]
		for item in data:
			gene = item.split("*")[0]
			allele = item
			if gene in spechla[platform] and sample in spechla[platform][gene]:
				spechla[platform][gene][sample].append(allele)

def populate_hla_resolve_dict():
	with open(hla_resolve, "r") as f:
		resolve_lines = f.read().splitlines()

	for line in resolve_lines[1:]:
		fields = line.split(",")
		sample = fields[0]
		tag = fields[1]  # keep full tag, e.g. "revio_hiphase"

		alleles = fields[2:]
		for allele in alleles:
			gene = allele.split("*")[0]
			if tag in resolve and gene in resolve[tag] and sample in resolve[tag][gene]:
				resolve[tag][gene][sample].append(allele)

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

def write_results():
	with open(output_file, "w", newline="") as f:
		writer = csv.writer(f)

		# Create header: sample, source, A_1, A_2, ..., DRB1_1, DRB1_2
		gene_pairs = [f"{gene}_1" for gene in genes] + [f"{gene}_2" for gene in genes]
		header = ["sample", "source"] + gene_pairs
		writer.writerow(header)

		for sample in samples:
			# IHW typing (ground truth)
			row = [sample, "IHW"]
			for gene in genes:
				alleles = ihw_data_dict.get(sample, {}).get("hla", {}).get(gene, ["", ""])
				alleles += [""] * (2 - len(alleles))
				if all("*" in a for a in alleles[:2]):
					alleles[:2] = sorted(alleles[:2], key=lambda x: int(x.split("*")[1].split(":")[0]))
				row.extend(alleles[:2])
			writer.writerow(row)

			# SpecHLA results
			for platform in spechla_platforms:
				row = [sample, f"{platform}_spechla"]
				for gene in genes:
					alleles = spechla[platform][gene].get(sample, [])
					alleles += [""] * (2 - len(alleles))
					if all("*" in a for a in alleles[:2]):
						alleles[:2] = sorted(alleles[:2], key=lambda x: int(x.split("*")[1].split(":")[0]))
					row.extend(alleles[:2])
				writer.writerow(row)

			# Resolve results from each phasing platform
			for tag in resolve_platforms:
				row = [sample, f"{tag}_resolve"]
				for gene in genes:
					alleles = resolve[tag][gene].get(sample, [])
					alleles += [""] * (2 - len(alleles))
					if all("*" in a for a in alleles[:2]):
						alleles[:2] = sorted(alleles[:2], key=lambda x: int(x.split("*")[1].split(":")[0]))
					row.extend(alleles[:2])
				writer.writerow(row)

def main():
	populate_spechla_dict()
	populate_hla_resolve_dict()
	populate_ihw_dict()
	write_results()

if __name__ == "__main__":
	main()
