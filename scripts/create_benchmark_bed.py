import os

# Original bed records for modification 
# chr6	31353871	31367067	HLA-B_ENSG00000234745
# chr6	29722774	29738528	HLA-F_ENSG00000204642
# chr6	32628178	32647062	HLA-DQA1_ENSG00000196735
# chr6	31982056	32002681	C4A_ENSG00000244731
# chr6	32014794	32035418	C4B_ENSG00000224389

# Regions that were not targeted and should be trimmed/masked

# HLA-F: [29728495,29738528]
# 29728495 is the BED stop of the last captured exon (ENSE00001697425). Cut everything downstream of this exon. 

# HLA-B: [31357681,31367067]
# 31357681 is BED stop of the first captured exon (ENSE00003967714; gene on negative strand). Cut everything upstream of of this exon. 

# C4A: [31984398,31991191]
# Remove intron with HERV insertion that causes reference mapping problems 

# C4B: [32017136,32023929]
# Remove intron with HERV insertion that causes reference mapping problems 

# HLA-DQA1: [32628178,32632743], [32643685,32647062]
# 32632743 is the BED start of the first captured exon (ENSE00001679862). Cut everything upstream.
# Cut everything downstream of the MANE 3' UTR BED stop (32643684)

root_dir = "/Users/matt/Documents/GitHub/mhc/bed_files/"

captured_genes_bed = root_dir + "hla_captured_genes.bed"
regions_to_exclude_bed = root_dir + "regions_to_exclude.bed"

output_dir = root_dir + "benchmark/"
make_output_dir = "mkdir -p {}".format(output_dir)
os.system(make_output_dir)

def strike_upcaptured_regions():
	subtract_regions_to_strike = "bedtools subtract -a {} -b {} > {}hla_captured_genes_modified.bed".format(captured_genes_bed, regions_to_exclude_bed, output_dir)
	merge_overlapping_intervals = "bedtools sort -i {}hla_captured_genes_modified.bed | bedtools merge -i - > {}merged_hla.bed".format(output_dir, output_dir)

	os.system(subtract_regions_to_strike)
	os.system(merge_overlapping_intervals)

def main():
	strike_upcaptured_regions()

if __name__ == "__main__":
	main()