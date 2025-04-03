import os 
import subprocess
import json

our_pacbio_seqs = "pacbio_fasta_dict.json"
our_ont_seqs = "ont_fasta_dict.json"

pacbio_hla_call_dir = "/hb/scratch/ogarci12/SPECHLA_HC/PacBio_HLA_Calls"
ont_hla_call_dir = "/hb/scratch/ogarci12/SPECHLA_HC/ONT_HLA_Calls"

samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", "IHW09409", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]
genes_of_interest = ["HLA-A", "HLA-B", "HLA-C"]

pacbio_seq_dict = {gene: {sample: [] for sample in samples} for gene in genes_of_interest}
ont_seq_dict = {gene: {sample: [] for sample in samples} for gene in genes_of_interest}

def get_fasta_files(dir):
	find_cmd = "find {input_dir} -type f -name *.fasta* | grep -v 'raw' | grep -v 'fai' > fasta_files.txt".format(input_dir = dir)
	subprocess.run(find_cmd, shell = True, check = True)
	with open("fasta_files.txt", "r") as f:
		files = f.read().splitlines()
	os.remove("fasta_files.txt")
	return files

def parse_fasta(platform_dict, file):
	sample = file.split("/")[-2]
	fields = file.split("/")[-1]
	allele_num = fields.split(".")[2]
	gene = fields.split(".")[3].replace("_", "-")
	if sample in samples and gene in genes_of_interest:

		with open(file, "r") as f:
			lines = f.read().split(">")

		seq = lines[1].split("\n",1)[1].strip().replace("\n", "")

		platform_dict[gene][sample].append(seq)

def combine_dicts(our_pacbio_seqs):
# def combine_dicts(our_pacbio_seqs, our_ont_seqs):
	with open(our_pacbio_seqs, "r") as f:
		our_pacbio = json.load(f)
	# with open(out_ont_seqs, "r") as f:
	# 	out_ont = json.load(f)

	for gene in our_pacbio:
		for sample in our_pacbio[gene]:
			output_file = f"{sample}_{gene}.fa"
			with open(output_file, "w") as f:
				f.write(f">{sample}_{gene}_pacbio_1" + "\n")
				f.write(our_pacbio[gene][sample][0] + "\n")
				
				f.write(f">{sample}_{gene}_pacbio_2" + "\n")
				f.write(our_pacbio[gene][sample][1] + "\n")

				# f.write(f">{sample}_{gene}_ont_1" + "\n")
				# f.write(our_ont[gene][sample][0] + "\n")
				
				# f.write(f">{sample}_{gene}_ont_2" + "\n")
				# f.write(our_ont[gene][sample][1] + "\n")
				
				f.write(f">{sample}_{gene}_pacbio_specHLA_1" + "\n")
				f.write(pacbio_seq_dict[gene][sample][0])
				
				f.write(f">{sample}_{gene}_pacbio_specHLA_2" + "\n")
				f.write(pacbio_seq_dict[gene][sample][1] + "\n")

				f.write(f">{sample}_{gene}_ont_specHLA_1" + "\n")
				f.write(ont_seq_dict[gene][sample][0] + "\n")
				
				f.write(f">{sample}_{gene}_ont_specHLA_2" + "\n")
				f.write(ont_seq_dict[gene][sample][1] + "\n")

def main():
	pacbio_fasta_files = get_fasta_files(pacbio_hla_call_dir)
	ont_fasta_files = get_fasta_files(ont_hla_call_dir)

	for file in pacbio_fasta_files:
		parse_fasta(pacbio_seq_dict, file)

	for file in ont_fasta_files:
		parse_fasta(ont_seq_dict, file)

	combine_dicts(our_pacbio_seqs)
	# combine_dicts(our_pacbio_seqs, our_ont_seqs)

if __name__ == "__main__":
	main()
