import os
import csv
import pysam

# HLA gene coordinates
genes_dict = {
	"HLA-A": [29941259, 29949572],
	"HLA-C": [31268748, 31272130],
	"HLA-B": [31353871, 31357681],
	"HLA-DRB5": [32517352, 32530287],
	"HLA-DRB1": [32577901, 32589848],
	"HLA-DQA1": [32632743, 32643685],
	"HLA-DQB1": [32659466, 32668383],
	"HLA-DQA2": [32741390, 32747198],
	"HLA-DQB2": [32756097, 32763532],
	"HLA-DPA1": [33064568, 33080775],
	"HLA-DPB1": [33075989, 33089696]
}

genes_of_interest = list(genes_dict.keys())

# Downsample directories (format: 0.005, 0.01, ..., 0.1)
downsample_dirs = [f"{i/1000:.3f}".rstrip("0").rstrip(".") for i in range(5, 105, 5)]

phased_results = []

for ds_dir in downsample_dirs:
	sample_results = [ds_dir]

	for gene in genes_of_interest:
		g_start, g_end = genes_dict[gene]
		g_len = g_end - g_start

		# Choose class based on gene name
		if gene in {"HLA-A", "HLA-B", "HLA-C"}:
			prefix = "MHC_Class_I"
		else:
			prefix = "MHC_Class_II"

		haploblock_file = os.path.join(ds_dir, f"{prefix}.phased.blocks.txt")
		vcf_file = os.path.join(ds_dir, f"{prefix}.dedup.trimmed.hg38.chr6.phased.vcf.gz")

		if not os.path.isfile(haploblock_file) or not os.path.isfile(vcf_file):
			print(f"Missing file(s) for {gene} in {ds_dir}")
			sample_results.append(0)
			continue

		# Load heterozygous variants
		het_sites = []
		vcf = pysam.VariantFile(vcf_file)
		for record in vcf:
			if record.chrom != "chr6":
				continue
			gt = record.samples["HG002"].get("GT")
			if gt in [(0,1), (1,0)]:
				het_sites.append(record.pos)

		gene_hets = [h for h in het_sites if g_start <= h <= g_end]
		print(f"{ds_dir} {gene}: {len(gene_hets)} hets, length={g_len}")

		if len(gene_hets) <= 1:
			sample_results.append(1)
			continue

		# Load haploblocks
		haploblocks = []
		with open(haploblock_file) as f:
			next(f)
			for line in f:
				fields = line.strip().split("\t")
				chrom, start, end = fields[3], int(fields[4]), int(fields[5])
				if chrom == "chr6":
					haploblocks.append((start, end))

		fully_phased = False
		for h_start, h_end in haploblocks:
			if h_start <= g_start and h_end >= g_end:
				fully_phased = True
				break
			if h_start <= g_end and h_end >= g_start:
				ext_start = max([h for h in het_sites if h < h_start], default=h_start)
				ext_start = max(ext_start, g_start)
				ext_end = min([h for h in het_sites if h > h_end], default=h_end)
				ext_end = min(ext_end, g_end)
				if ext_start <= g_start and ext_end >= g_end:
					fully_phased = True
					break

		sample_results.append(1 if fully_phased else 0)

	phased_results.append(sample_results)

# Write CSV
output_file = "hla_class_i_ii_phasing_summary.csv"
with open(output_file, "w", newline="") as f:
	writer = csv.writer(f)
	writer.writerow(["sample"] + genes_of_interest)
	writer.writerows(phased_results)

print(f"Summary written to: {output_file}")
