import pysam
import csv

# Input paths
vcf_path = "/hb/groups/cornejo_lab/matt/hla_capture/pacbio/phased_vcf_hiphase/HG002.dedup.trimmed.hg38.chr6.phased.joint.vcf.gz"
bed_path = "mhc_genes.bed"
output_counts = "unphased_het_counts.tsv"

# Load MHC regions and preserve original order
gene_regions = []
gene_order = []
with open(bed_path) as f:
    for line in f:
        chrom, start, end, full_gene = line.strip().split("\t")
        gene_name = full_gene.split("_")[0]  # Strip ENSG part
        gene_regions.append((chrom, int(start), int(end), gene_name))
        if gene_name not in gene_order:
            gene_order.append(gene_name)

# Initialize counts in original gene order
gene_counts = {gene: 0 for gene in gene_order}

# Open VCF
vcf = pysam.VariantFile(vcf_path)
region = "chr6"
start = 28000000
end = 34000000

# Count unphased hets per gene
for record in vcf.fetch(region, start, end):
    for sample in record.samples:
        gt = record.samples[sample]['GT']
        phased = record.samples[sample].phased
        if gt and len(gt) == 2 and gt[0] != gt[1] and not phased:
            pos = record.pos - 1  # BED-style
            for chrom, r_start, r_end, gene in gene_regions:
                if chrom == record.chrom and r_start <= pos < r_end:
                    gene_counts[gene] += 1
                    break  # only count once per variant

# Write output in original BED order
with open(output_counts, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["Gene", "Unphased_Het_Count"])
    for gene in gene_order:
        writer.writerow([gene, gene_counts[gene]])
