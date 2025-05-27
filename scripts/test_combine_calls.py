import os
import csv

# Core setup
genes = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
samples = [
	"HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055",
	"HG02630", "HG03492", "HG03579", "IHW09021", "IHW09049", "IHW09071", "IHW09117",
	"IHW09118", "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", "IHW09224",
	"IHW09245", "IHW09251", "IHW09359", "IHW09364", "IHW09409", "NA19240", "NA20129",
	"NA21309", "NA24694", "NA24695"
]

# File and directory paths
spechla_pacbio = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/spechla_results/PacBio_spechla_results.txt"
spechla_ont = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/spechla_results/ONT_spechla_results.txt"
hla_resolve = "hla_typing_results.csv"
ihw_data_file = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/IPD_IMGT/IHW_samples.txt"
specimmune_dir = "/hb/scratch/ogarci12/SpecImmune/"
output_file = "typing_comparison.csv"

# Platforms
spechla_platforms = ["revio", "promethion"]
resolve_platforms = ["revio_hiphase", "revio_longphase", "promethion_longphase"]
specimmune_platforms = ["PacBio_HLA_Calls", "ONT_HLA_Calls"]

# Data containers
spechla = {p: {g: {s: [] for s in samples} for g in genes} for p in spechla_platforms}
resolve = {p: {g: {s: [] for s in samples} for g in genes} for p in resolve_platforms}
specimmune = {p.lower(): {g: {s: [] for s in samples} for g in genes} for p in specimmune_platforms}
ihw_data_dict = {}

hla_la_dir = "/hb/scratch/ogarci12/HLA_HC_Call/HLA-LA/"
hla_la_platforms = ["PB", "ONT"]
hla_la = {p.lower(): {g: {s: [] for s in samples} for g in genes} for p in hla_la_platforms}


def populate_specimmune_dict():
	for platform in specimmune_platforms:
		platform_path = os.path.join(specimmune_dir, platform)
		tag = platform.lower()
		for sample in samples:
			result_file = os.path.join(platform_path, sample, f"{sample}.HLA.final.type.result.formatted.txt")
			if not os.path.isfile(result_file):
				continue
			with open(result_file, "r") as f:
				for line in f:
					if line.startswith("#") or not line.strip():
						continue
					fields = line.strip().split("\t")
					if len(fields) < 7:
						continue
					locus = fields[0].strip()
					one_guess = fields[6].strip().replace("HLA-", "")
					gene = locus.replace("HLA-", "")
					if gene in genes:
						specimmune[tag][gene][sample].append(one_guess)

def populate_spechla_dict():
	for path, platform in [(spechla_pacbio, "revio"), (spechla_ont, "promethion")]:
		with open(path, "r") as f:
			for line in f.readlines()[2:]:
				parts = line.strip().split("\t")
				sample, alleles = parts[0], parts[1:]
				for allele in alleles:
					gene = allele.split("*")[0]
					if gene in spechla[platform] and sample in spechla[platform][gene]:
						spechla[platform][gene][sample].append(allele)

def populate_hlala_dict():
	for platform in hla_la_platforms:
		base_path = os.path.join(hla_la_dir, platform)
		for sample in samples:
			bestguess_file = os.path.join(base_path, sample, sample, "hla", "R1_bestguess_G.txt")
			if not os.path.isfile(bestguess_file):
				continue
			with open(bestguess_file, "r") as f:
				next(f)  # skip header
				for line in f:
					fields = line.strip().split("\t")
					if len(fields) < 3:
						continue
					locus = fields[0].strip()
					allele = fields[2].strip()
					gene = locus.replace("HLA-", "")  # just in case
					if gene in genes:
						hla_la[platform.lower()][gene][sample].append(allele)

def populate_hla_resolve_dict():
	with open(hla_resolve, "r") as f:
		next(f)
		for line in f:
			sample, tag, *alleles = line.strip().split(",")
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

def normalize_alleles(alleles):
	alleles = [a.strip() for a in alleles if a.strip()]
	alleles += [""] * (2 - len(alleles))
	if all("*" in a for a in alleles[:2]):
		try:
			alleles[:2] = sorted(alleles[:2], key=lambda x: int(x.split("*")[1].split(":")[0]))
		except:
			pass
	return alleles[:2]

def write_results():
	with open(output_file, "w", newline="") as f:
		writer = csv.writer(f)

		# Header: sample, source, A_1, A_2, B_1, B_2, ...
		header = ["sample", "source"]
		for gene in genes:
			header.extend([f"{gene}_1", f"{gene}_2"])
		writer.writerow(header)

		for sample in samples:
			def write_row(source, data_dict):
				row = [sample, source]
				for gene in genes:
					alleles = normalize_alleles(data_dict[gene][sample])
					row.extend(alleles)
				writer.writerow(row)

			# IHW
			ihw_row = [sample, "IHW"]
			for gene in genes:
				alleles = normalize_alleles(ihw_data_dict.get(sample, {}).get("hla", {}).get(gene, []))
				ihw_row.extend(alleles)
			writer.writerow(ihw_row)

			# SpecHLA
			for platform in spechla_platforms:
				write_row(f"{platform}_spechla", spechla[platform])

			# Resolve
			for platform in resolve_platforms:
				write_row(f"{platform}_resolve", resolve[platform])

			# SpecImmune with renamed tags
			platform_map = {
				"pacbio_hla_calls": "revio_specimmune",
				"ont_hla_calls": "promethion_specimmune"
			}
			for platform in specimmune_platforms:
				tag = platform.lower()
				label = platform_map.get(tag, f"{tag}_specimmune")
				write_row(label, specimmune[tag])

			# HLA*LA
			hla_la_map = {
				"pb": "revio_hlala",
				"ont": "promethion_hlala"
			}
			for platform in hla_la_platforms:
				tag = platform.lower()
				label = hla_la_map.get(tag, f"{tag}_hlala")
				write_row(label, hla_la[tag])

def main():
	populate_spechla_dict()
	populate_hla_resolve_dict()
	populate_specimmune_dict()
	populate_hlala_dict()
	populate_ihw_dict()
	write_results()

if __name__ == "__main__":
	main()
