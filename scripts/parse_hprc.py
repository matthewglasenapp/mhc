import os
import gzip
import csv
import subprocess

samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]
genes = ["HLA-A", "HLA-B", "HLA-C"]

hprc_ensembl_fasta_dir = "/hb/scratch/mglasena/hprc_ensembl/"
hprc_ensembl_gff3_dir = "/hb/scratch/mglasena/hprc_gff3/"
# output_bed_dir = "/hb/scratch/mglasena/hprc_hla_beds/"
output_bed_dir = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/hprc/"
os.makedirs(output_bed_dir, exist_ok=True)

sample_to_gca = {
	'HG01258.pri.mat.f1_v2': 'GCA_018469405.1',
	'HG01258.alt.pat.f1_v2': 'GCA_018469675.1',
	'HG02630.alt.pat.f1_v2': 'GCA_018469945.1',
	'HG02630.pri.mat.f1_v2': 'GCA_018469955.1',
	'HG01106.alt.pat.f1_v2': 'GCA_018471075.1',
	'HG01106.pri.mat.f1_v2': 'GCA_018471345.1',
	'HG01928.pri.mat.f1_v2': 'GCA_018472695.1',
	'HG01928.alt.pat.f1_v2': 'GCA_018472705.1',
	'HG03579.pri.mat.f1_v2': 'GCA_018472825.1',
	'HG03579.alt.pat.f1_v2': 'GCA_018472835.1',
	'NA19240.alt.pat.f1_v2': 'GCA_018503265.1',
	'NA19240.pri.mat.f1_v2': 'GCA_018503275.1',
	'NA20129.alt.pat.f1_v2': 'GCA_018504625.1',
	'NA20129.pri.mat.f1_v2': 'GCA_018504635.1',
	'NA21309.pri.mat.f1_v2': 'GCA_018504655.1',
	'NA21309.alt.pat.f1_v2': 'GCA_018504665.1',
	'HG03492.alt.pat.f1_v2': 'GCA_018505835.1',
	'HG03492.pri.mat.f1_v2': 'GCA_018505845.1',
	'HG02055.alt.pat.f1_v2': 'GCA_018505855.1',
	'HG02055.pri.mat.f1_v2': 'GCA_018506125.1',
	'HG005.alt.pat.f1_v2': 'GCA_018506945.1',
	'HG005.pri.mat.f1_v2': 'GCA_018506965.1',
	# 'HG002.alt.pat.f1_v2': 'GCA_018852605.1',
	# 'HG002.pri.mat.f1_v2': 'GCA_018852615.1',
	'HG002.pat.cur.20211005': 'GCA_021950905.1',
	'HG002.mat.cur.20211005': 'GCA_021951015.1'
}

def extract_hla_coordinates(assembly_name, GCA_ID):
	gff3_file = os.path.join(hprc_ensembl_gff3_dir, f"{GCA_ID}-genes.gff3.gz")
	bed_file = os.path.join(output_bed_dir, f"{assembly_name}_hla.bed")

	with gzip.open(gff3_file, "rt") as infile, open(bed_file, "w", newline="") as outbed:
		writer = csv.writer(outbed, delimiter="\t")
		for line in infile:
			if line.startswith("#"):
				continue
			fields = line.strip().split("\t")
			
			if fields[2] == "gene":
				info = fields[8]
				
				for gene in genes:
					if f"Name={gene}" in info:
						chrom = fields[0]
						# Convert to BED
						start = str(int(fields[3]) - 1)
						end = fields[4]
						strand = fields[6]
						writer.writerow([chrom, start, end, gene, ".", strand])
						break
	
	return bed_file

def extract_fasta(assembly_name, GCA_ID, bed_file):
	fasta_file = os.path.join(hprc_ensembl_fasta_dir, f"Homo_sapiens-{GCA_ID}-unmasked.fa.gz")
	output_fasta = os.path.join(output_bed_dir, f"{assembly_name}_hla.fa")
	cmd = ["bedtools", "getfasta", "-fi", fasta_file, "-bed", bed_file, "-fo", output_fasta, "-name", "-s"]
	subprocess.run(cmd, check=True)

def main():
	for sample in samples:
		for assembly_name, GCA_ID in sample_to_gca.items():
			if sample == assembly_name.split(".")[0]:
				print(f"🧬 Processing {assembly_name}")
				bed_file = extract_hla_coordinates(assembly_name, GCA_ID)
				extract_fasta(assembly_name, GCA_ID, bed_file)

if __name__ == "__main__":
	main()

