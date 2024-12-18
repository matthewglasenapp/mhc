import os
from joblib import Parallel, delayed

bam_file_dir = "/hb/scratch/mglasena/MHC/mapped_bam/"
output_dir = "/hb/scratch/mglasena/MHC/coverage/"

regions_file = "/hb/scratch/mglasena/MHC/scripts/mhc_regions.bed"

threads = 4
mapq_threshold = 20

def get_bam_file_paths():
	get_bam_paths_file = 'find {} -type f -name *.bam* | grep -v "bai" > bam_files.txt'.format(bam_file_dir)
	os.system(get_bam_paths_file)

	with open("bam_files.txt", "r") as f:
		bam_file_paths_list = f.read().splitlines()

	os.system("rm bam_files.txt")

	return bam_file_paths_list

def index_bam(bam_file):
	index = "samtools index -b {}".format(bam_file)
	os.system(index)

def run_mosdepth(bam_file):
	os.chdir(output_dir)
	prefix = bam_file.split("/")[-1].split(".")[0] + "_" + bam_file.split(".")[2].split("_")[0]
	# --flag 3328 excludes duplicates and secondary/supplementary alignments
	mosdepth = "mosdepth --flag 3328 --mapq {} --by {} --thresholds 20,30 -t {} {} {}".format(mapq_threshold, regions_file, threads, prefix, bam_file)
	os.system(mosdepth)

def main():
	bam_file_paths_list = get_bam_file_paths()
	#Parallel(n_jobs=24)(delayed(index_bam)(bam_file) for bam_file in bam_file_paths_list)
	Parallel(n_jobs=6)(delayed(run_mosdepth)(bam_file) for bam_file in bam_file_paths_list)

if __name__ == "__main__":
	main()