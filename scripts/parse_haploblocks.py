import os
import csv
import json
from joblib import Parallel, delayed

haploblock_dir = "/hb/scratch/mglasena/MHC/scripts/haploblocks/"
genes_bed = "hla_captured_genes.bed"
output_csv = "phased_genes.tsv"
output_json = "gene_haploblock_dict.json"

samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01891", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", "IHW09409", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]

genes_dict = dict()
haploblock_dict = {sample: [] for sample in samples}
gene_haploblock_dict = {sample: [] for sample in samples}

hla_start = 29722774
hla_stop = 33129084

def create_genes_dict():
    genes = open(genes_bed, "r").read().splitlines()
    for line in genes:
        fields = line.split("\t")
        name = fields[3].split("_")[0]
        start = int(fields[1])
        stop = int(fields[2])
        genes_dict[name] = [start, stop]

def parse_haploblocks(sample):
    platforms = ["revio", "promethion"]
    haploblock_lists = []

    for platform in platforms:
        haploblock_file = haploblock_dir + sample + "_" + platform + "_haploblocks.gtf"
        haploblocks = open(haploblock_file, "r").read().splitlines()
        platform_haploblock_list = []
        for line in haploblocks:
            fields = line.split("\t")
            start = int(fields[3]) - 1
            stop = int(fields[4])
            if stop > hla_start:
                platform_haploblock_list.append([start, stop])
        haploblock_lists.append(platform_haploblock_list)

    return sample, haploblock_lists[0], haploblock_lists[1]

def evaluate_gene_haploblocks(sample):
    revio_gene_list = []
    promethion_gene_list = []
    revio_haploblocks = haploblock_dict[sample][0]
    promethion_haploblocks = haploblock_dict[sample][1]

    #print(f"\n[Revio Platform Results for {sample}]:")
    for gene in genes_dict:
        gene_start = genes_dict[gene][0]
        gene_stop = genes_dict[gene][1]

        # Check Revio blocks
        for block_start, block_stop in revio_haploblocks:
            if block_start <= gene_start and block_stop >= gene_stop:
                #print(f"Block: {block_start}-{block_stop}, Gene: {gene_start}-{gene_stop}")
                #print(f"Gene {gene} is fully contained in block {block_start}-{block_stop}")
                revio_gene_list.append(gene)
                break

    #print(f"\n[PromethION Platform Results for {sample}]:")
    for gene in genes_dict:
        gene_start = genes_dict[gene][0]
        gene_stop = genes_dict[gene][1]

        # Check PromethION blocks
        for block_start, block_stop in promethion_haploblocks:
            if block_start <= gene_start and block_stop >= gene_stop:
                #print(f"Block: {block_start}-{block_stop}, Gene: {gene_start}-{gene_stop}")
                #print(f"Gene {gene} is fully contained in block {block_start}-{block_stop}")
                promethion_gene_list.append(gene)
                break

    return sample, revio_gene_list, promethion_gene_list

def write_results():
    with open(output_csv, "w", newline="") as csv_file:
        writer = csv.writer(csv_file, delimiter="\t")
        writer.writerow(["sample", "platform", "num_genes", "genes"])
        for sample, gene_lists in gene_haploblock_dict.items():
            revio_genes = gene_lists[0]
            promethion_genes = gene_lists[1]

            writer.writerow([sample, "Revio", len(revio_genes), ",".join(revio_genes)])
            writer.writerow([sample, "PromethION", len(promethion_genes), ",".join(promethion_genes)])

    platform_dict = {sample: {"Revio": gene_lists[0], "PromethION": gene_lists[1]} for sample, gene_lists in gene_haploblock_dict.items()}
    with open(output_json, "w") as json_file:
        json.dump(platform_dict, json_file, indent=4)

def main():
    create_genes_dict()

    haploblocks_by_sample = Parallel(n_jobs=10)(delayed(parse_haploblocks)(sample) for sample in samples)

    for sample, revio_haploblock_list, promethion_haploblock_list in haploblocks_by_sample:
        haploblock_dict[sample] = [revio_haploblock_list, promethion_haploblock_list]

    genes_by_haploblock = Parallel(n_jobs=10)(delayed(evaluate_gene_haploblocks)(sample) for sample in samples)

    for sample, revio_gene_list, promethion_gene_list in genes_by_haploblock:
        gene_haploblock_dict[sample] = [revio_gene_list, promethion_gene_list]

    write_results()

if __name__ == "__main__":
    main()

