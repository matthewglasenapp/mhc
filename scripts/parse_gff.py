import os
import gzip
import json

root_dir = "/Users/matt/Documents/GitHub/mhc/"

output_dir = root_dir + "bed_files/"
make_output_dir = "mkdir -p {}".format(output_dir)
os.system(make_output_dir)

raw_data_dir = root_dir + "raw_data/"
make_raw_data_dir = "mkdir -p {}".format(raw_data_dir)
os.system(make_raw_data_dir)

hg38_gff_ftp = "https://ftp.ensembl.org/pub/release-110/gff3/homo_sapiens/Homo_sapiens.GRCh38.110.chromosome.6.gff3.gz"

# Initialize dictionary of {transcript_id: parent_gene_name}
transcript_parent_gene_dict = dict()

# Gene record: biotype
biotype_dict = dict()

# Record types to include in trasncript_parent_gene_dict
record_types = ['lnc_rna', 'mirna', 'mrna', 'ncrna', 'pseudogenic_transcript', 'snorna', 'snrna', 'transcript', 'unconfirmed_transcript']

mhc_start = 28000000
mhc_stop = 34000000

def download_gff():
	# Get human GFF3 file for Genome Reference Consortium human genome build 38, release 110
	os.system(f"wget {hg38_gff_ftp} -P {raw_data_dir}")

def build_transcript_parent_gene_dict(gff_file):
	seen_transcripts = set()
	lines = gzip.open(gff_file,"rt").read().splitlines()
	for line in lines:
		
		# Ignore info lines
		if line[0] == "#":
			continue
		
		# Parse gff mRNA records for transcript ID and parent gene name
		record_type = line.split("\t")[2].lower()
		if record_type in record_types:
			info_field = line.split("\t")[8]

			transcript_id = None
			gene_name = None
			
			if "ID=transcript:" in info_field:
				transcript_id = info_field.split("ID=transcript:")[1].split(";")[0]
			if "Parent=gene:" in info_field:
				gene_id = info_field.split("Parent=gene:")[1].split(";")[0]
			if "Name=" in info_field:
				gene_name = info_field.split("Name=")[1].split(";")[0]
			else:
				print("No name")
				print(info_field)
				gene_name = gene_id

			if transcript_id:
			    if transcript_id in seen_transcripts:
			        print(f"Duplicate detected in input for {transcript_id}.")
			    seen_transcripts.add(transcript_id)
			
			if transcript_id in transcript_parent_gene_dict:
				print(f"Transcript {transcript_id} already in dict.")
				print(f"Existing entry: {transcript_parent_gene_dict[transcript_id]}")
			elif transcript_id and gene_name:
				#print(f"Adding {transcript_id} -> {gene_name} to transcript_parent_gene_dict.")
				transcript_parent_gene_dict[transcript_id] = gene_name
			else:
				print(f"Skipped line with missing data: {info_field}")

def create_bed_files(gff_file, output_file_genes, output_file_exons):
	genes_bed_record_types = ["gene", "ncrna_gene", "pseudogene"]
	gff_lines = gzip.open(gff_file,"rt").read().splitlines()

	with open(output_dir + output_file_genes, "w") as genes_bed, open(output_dir + output_file_exons, "w") as exons_bed:
		for line in gff_lines:

			# Ignore info lines
			if line[0] == "#":
				continue

			chromosome = line.split("\t")[0]
			if chromosome == "6":
				start = int(line.split("\t")[3]) - 1
				stop = int(line.split("\t")[4])
				if start >= mhc_start and stop <= mhc_stop:
					record_type = line.split("\t")[2].lower()
					if record_type == "exon" or record_type in genes_bed_record_types:
						chromosome = "chr" + chromosome
						start = str(start)
						stop = str(stop)
						info_field = line.split("\t")[8]
					
					if record_type in genes_bed_record_types:
						
						gene_id = info_field.split("ID=gene:")[1].split(";")[0]
						biotype = info_field.split("biotype=")[1].split(";")[0]
						biotype_dict[gene_id] = biotype
						
						if "Name=" in info_field:
							gene_name = info_field.split("Name=")[1].split(";")[0]
							name_string = gene_name + "_" + gene_id
						
						else:
							name_string = gene_id

						genes_bed.write(chromosome + "\t" + start + "\t" + stop + "\t" + name_string + "\n")

					elif record_type == "exon":
						exon_id = info_field.split("exon_id=")[1].split(";")[0]
						transcript_id = info_field.split("Parent=transcript:")[1].split(";")[0].strip()
						parent_gene = transcript_parent_gene_dict[transcript_id]
						name_string = exon_id + "_" + transcript_id + "_" + parent_gene
						exons_bed.write(chromosome + "\t" + start + "\t" + stop + "\t" + name_string + "\n")

	with open(output_dir + "biotype_dict.json", "w") as json_file:
		json.dump(biotype_dict, json_file, indent = 4)

def create_mosdepth_regions_file(output_file):
	with open(output_dir + "mhc_genes.bed", "r") as f1, open(output_dir + "mhc_exons.bed", "r") as f2:
		genes = f1.read().splitlines()
		exons = f2.read().splitlines()

	with open(output_dir + output_file, "w") as f3:
		for line in genes:
			fields = line.split("\t")
			new_field_four = ["gene_" + fields[3]]
			new_line = "\t".join(fields[0:3] + new_field_four)
			f3.write(new_line + "\n")
		
		for line in exons:
			fields = line.split("\t")
			new_field_four = ["exon_" + fields[3]]
			new_line = "\t".join(fields[0:3] + new_field_four)
			f3.write(new_line + "\n")
					
def main():
	#download_gff()
	gff3_file = "{}Homo_sapiens.GRCh38.110.chromosome.6.gff3.gz".format(raw_data_dir)
	build_transcript_parent_gene_dict(gff3_file)
	create_bed_files(gff3_file, "mhc_genes.bed", "mhc_exons.bed")
	create_mosdepth_regions_file("mhc_regions.bed")

if __name__ == "__main__":
	main()