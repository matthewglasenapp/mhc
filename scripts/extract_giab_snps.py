import os
import csv

# Input files
vcf_file = "/hb/scratch/mglasena/MHC/concordance/GIAB_benchmark/HG002/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
confident_regions_bed = "/hb/scratch/mglasena/MHC/concordance/GIAB_benchmark/HG002/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
captured_mhc_genes_bed = "hla_captured_genes.bed"

output_dir = "/hb/scratch/mglasena/MHC/concordance/GIAB_SNP/"
make_output_dir = "mkdir -p {}".format(output_dir)
os.system(make_output_dir)

sample = "HG002"

# Output file
output_csv = output_dir + sample + "_gene_snps.csv"

def extract_snps():
    # 1. Extract SNPs from the VCF file.
    print("Extracting SNPs from {}".format(vcf_file))
    output_snp = output_dir + sample + "_SNP.vcf.gz"
    extract_snps = "bcftools view -v snps {} -Oz -o {}".format(vcf_file, output_snp)
    index = "bcftools index {}".format(output_snp)
    os.system(extract_snps)
    os.system(index)

def intersect_with_confident_regions():
    # 2. Intersect SNPs with confident regions.
    snp_file = output_dir + sample + "_SNP.vcf.gz"
    output_file = output_dir + sample + "_SNP_confident.vcf.gz"
    print("Intersecting {} with {}".format(snp_file, confident_regions_bed))
    intersect = "bedtools intersect -a {} -b {} -header | bgzip > {}".format(snp_file, confident_regions_bed, output_file)
    index = "bcftools index {}".format(output_file)
    os.system(intersect)
    os.system(index)

def intersect_with_genes():
    # 3. Intersect confident regions SNPs with MHC genes.
    confident_snp_file = output_dir + sample + "_SNP_confident.vcf.gz"
    output_file = output_dir + sample + "_SNP_confident_genes.vcf.gz"
    print("Intersecting {} with {}".format(confident_snp_file, captured_mhc_genes_bed))
    intersect = "bedtools intersect -a {} -b {} -header | bgzip > {}".format(confident_snp_file, captured_mhc_genes_bed, output_file)
    os.system(intersect)

def count_snps_per_gene():
    # 4. Count SNPs per gene.
    print("Calculating SNPs per gene")
    confident_snp_genes_file = output_dir + sample + "_SNP_confident_genes.vcf.gz"
    output_file = output_dir + sample + "_gene_snps.bed"
    intersect = "bedtools intersect -a {} -b {} -c > {}".format(captured_mhc_genes_bed, confident_snp_genes_file, output_file)
    os.system(intersect)

def convert_bed_to_csv():
    # 5. Convert the BED file with SNP counts to a CSV file.
    
    print("Converting BED file to CSV format")
    bed_file = output_dir + sample + "_gene_snps.bed"
    with open(bed_file, "r") as bed_file, open(output_csv, "w", newline="") as csv_file:
        csv_writer = csv.writer(csv_file)
        # Write the header row
        csv_writer.writerow(["Chromosome", "Start", "End", "Gene", "SNPs"])
        for line in bed_file:
            fields = line.strip().split("\t")
            chromosome = fields[0]
            start = fields[1]
            end = fields[2]
            gene_name = fields[3]
            snp_count = fields[-1]
            csv_writer.writerow([chromosome, start, end, gene_name, snp_count])

def main():    
    extract_snps()
    intersect_with_confident_regions()
    intersect_with_genes()
    count_snps_per_gene()
    convert_bed_to_csv()
    
if __name__ == "__main__":
    main()
