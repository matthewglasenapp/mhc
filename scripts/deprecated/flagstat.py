import os
from joblib import Parallel, delayed

bam_file_dir = "/hb/scratch/mglasena/MHC/mapped_bam/"

def get_file_paths_list():
	get_bam_paths_file = 'find {} -type f -name *.bam* | grep -v "bai" | grep -v "rmdup" > bam_files.txt'.format(bam_file_dir)
	os.system(get_bam_paths_file)

	with open("bam_files.txt", "r") as f:
		bam_file_paths_list = f.read().splitlines()

	#os.system("rm bam_files.txt")

	return bam_file_paths_list

def run_flagstat(bam_file):
	index = "samtools index {}".format(bam_file)
	os.system(index)
	output_file = bam_file.split(".")[0] + "_flagstat.tsv"
	flagstat = "samtools view -F 2304 -u {} | samtools flagstat -O tsv - > {}".format(bam_file, output_file)
	os.system(flagstat)

def main():
	bam_file_paths_list = get_file_paths_list()

	Parallel(n_jobs=len(bam_file_paths_list))(delayed(run_flagstat)(bam_file) for bam_file in bam_file_paths_list)

if __name__ == "__main__":
	main()
	