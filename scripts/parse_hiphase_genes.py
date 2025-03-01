import csv

phased_genes_file = "phased_genes.tsv"

samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", "IHW09409", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]
genes_of_interest = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2")

def parse_phased_genes():
	HLA_phased_dict = {gene: 0 for gene in genes_of_interest}

	with open(phased_genes_file,"r") as f:
		records = f.read().splitlines()[1:]

	for item in records:
		fields = item.split("\t")
		gene_list = fields[2].split(",")
		for gene in genes_of_interest:
			if gene in gene_list:
				HLA_phased_dict[gene] += 1

	return HLA_phased_dict

def make_heat_map_data():
	phased_dict = {sample: [] for sample in samples}
	
	with open(phased_genes_file,"r") as f:
		records = f.read().splitlines()[1:]

	for item in records:
		fields = item.split("\t")
		sample = fields[0]
		gene_list = fields[2].split(",")
		for gene in genes_of_interest:
			if gene in gene_list:
				phased_dict[sample].append(1)
			else:
				phased_dict[sample].append(0)

	return phased_dict

def write_results(phase_dict):
    output_csv = "phase_map.csv"
    with open(output_csv, "w", newline="") as csv_file:
        writer = csv.writer(csv_file, delimiter=",")
        header = ["sample", "HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2"]
        writer.writerow(header)
        for sample, genes in phase_dict.items():
        	writer.writerow([sample] + genes)

def main():
	HLA_phased_dict = parse_phased_genes()
	#print(HLA_phased_dict)
	prop_phased_dict = {gene: round(count / len(samples), 2) for gene, count in HLA_phased_dict.items()}
	print(prop_phased_dict)

	heat_map_data = make_heat_map_data()
	write_results(heat_map_data)

if __name__ == "__main__":
	main()