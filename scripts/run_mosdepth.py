import os
from joblib import Parallel, delayed

bam_file_dir = "/hb/scratch/ogarci12/deepvariant/"

regions_file = "mhc_regions.bed"

threads = 4

def get_bam_file_paths():
	get_files = 'find {} -type f -name "*.bam*" | grep -v "bai" > bam_file_paths.txt'.format(bam_file_dir)
	os.system(get_files)

	with open("bam_file_paths.txt","r") as f:
		bam_file_paths_list = f.read().splitlines()

	remove_bam_file_paths_file = "rm bam_file_paths.txt"
	os.system(remove_bam_file_paths_file)
	return bam_file_paths_list

def index_bam(bam_file):
	index = "samtools index -b {}".format(bam_file)
	os.system(index)

def run_mosdepth(bam_file):
	prefix = bam_file.split("/")[-1].split(".")[0] + "_" + bam_file.split("/")[-2].lower()
	mosdepth = "mosdepth --by {} --thresholds 20,30 -t {} {} {}".format(regions_file, threads, prefix, bam_file)
	os.system(mosdepth)

def main():
	#bam_file_paths_list = get_bam_file_paths()
	#Parallel(n_jobs=24)(delayed(index_bam)(bam_file) for bam_file in bam_file_paths_list)
	#Parallel(n_jobs=6)(delayed(run_mosdepth)(bam_file) for bam_file in bam_file_paths_list)

if __name__ == "__main__":
	main()