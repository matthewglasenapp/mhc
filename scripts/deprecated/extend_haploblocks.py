import pysam

# Input files
blocks_file = "HG002.phased.blocks.txt"
vcf_file = "HG002.dedup.trimmed.hg38.chr6.phased.vcf.gz"
output_bed = "HG002.extended.blocks.bed"

# Load heterozygous variants from HiPhase VCF
def load_heterozygous_variants():
	heterozygous_sites = {"chr6": []}
	vcf = pysam.VariantFile(vcf_file)
	for record in vcf:
		if record.chrom != "chr6":
			continue

		sample = list(record.samples.keys())[0]
		genotype = record.samples[sample]["GT"]

		if genotype in [(0, 1), (1, 0)]:
			heterozygous_sites["chr6"].append(record.pos)
	
	return heterozygous_sites 

# Extend phase blocks
def extend_phase_blocks(heterozygous_sites):
	extended_blocks = []

	with open(blocks_file) as f:
		lines = f.read().splitlines()
		for line in lines[1:]:
			fields = line.strip().split("\t")
			chrom = fields[3]
			start = int(fields[4])
			end = int(fields[5])

			if chrom != "chr6":
				continue
			
			prev_het = max([h for h in heterozygous_sites["chr6"] if h < start], default=start)
			next_het = min([h for h in heterozygous_sites["chr6"] if h > end], default=end)

			# Extend and convert to 0-based
			extended_blocks.append((chrom, prev_het - 1, next_het))
	
	return extended_blocks

# Write to BED file
def write_bed_file(extended_blocks):
	with open(output_bed, "w") as out:
		for chrom, start, end in extended_blocks:
			out.write(f"{chrom}\t{start}\t{end}\n")
	print(f"Extended phase blocks written to {output_bed}")

# Main function
def main():
	heterozygous_sites = load_heterozygous_variants()
	extended_blocks = extend_phase_blocks(heterozygous_sites)
	write_bed_file(extended_blocks)

if __name__ == "__main__":
	main()
