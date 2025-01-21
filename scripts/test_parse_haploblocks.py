import os
import csv
import json
from joblib import Parallel, delayed

genes_bed = "hla_captured_genes.bed"
output_csv = "phased_genes.csv"
output_json = "gene_haploblock_dict.json"

samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01891", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", "IHW09409", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]

genes_dict = dict()
haploblock_dict = {sample: [] for sample in samples}
gene_haploblock_dict = {sample: [] for sample in samples}

hla_start = 29722774
hla_stop = 33129084

def parse_haploblocks(sample, platform):
    haploblock_file = f"{sample}_{platform}_haploblocks.tsv"
    if not os.path.exists(haploblock_file):
        print(f"File not found: {haploblock_file}")
        return []

    haploblocks = open(haploblock_file, "r").read().splitlines()
    haploblock_list = []
    for line in haploblocks:
        start = int(line.split("\t")[3]) - 1
        stop = int(line.split("\t")[4])

        if stop > hla_start:
            haploblock_list.append([start, stop, platform])

    return haploblock_list

def create_genes_dict():
    genes = open(genes_bed, "r").read().splitlines()
    for line in genes:
        fields = line.split("\t")
        name = fields[3].split("_")[0]
        start = int(fields[1])
        stop = int(fields[2])
        genes_dict[name] = [start, stop]

def evaluate_gene_haploblocks(sample):
    gene_list = []
    haploblocks = haploblock_dict[sample]
    for gene in genes_dict:
        gene_start = genes_dict[gene][0]
        gene_stop = genes_dict[gene][1]

        for block_start, block_stop, platform in haploblocks:
            if block_start <= gene_start and block_stop >= gene_stop:
                gene_list.append((gene, platform))
                break

    return sample, gene_list

def write_results():
    with open(output_csv, "w", newline="") as csv_file:
        writer = csv.writer(csv_file, delimiter="\t")
        writer.writerow(["sample", "platform", "num_genes", "genes"])
        for sample, genes in gene_haploblock_dict.items():
            platform_dict = {}
            for gene, platform in genes:
                platform_dict.setdefault(platform, []).append(gene)

            for platform, gene_list in platform_dict.items():
                writer.writerow([sample, platform, len(gene_list), ",".join(gene_list)])

    with open(output_json, "w") as json_file:
        json.dump(gene_haploblock_dict, json_file, indent=4)

def main():
    haploblocks_by_sample = Parallel(n_jobs=len(samples))(
        delayed(parse_haploblocks)(sample, platform)
        for sample in samples
        for platform in ["revio", "promethion"]
    )

    for sample, platform_haploblocks in zip(samples * 2, haploblocks_by_sample):
        haploblock_dict[sample].extend(platform_haploblocks)

    create_genes_dict()

    genes_by_haploblock = Parallel(n_jobs=len(samples))(delayed(evaluate_gene_haploblocks)(sample) for sample in samples)

    for sample, gene_list in genes_by_haploblock:
        gene_haploblock_dict[sample] = gene_list

    write_results()

if __name__ == "__main__":
    main()
