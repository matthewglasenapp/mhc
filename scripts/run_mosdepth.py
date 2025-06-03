import os
from joblib import Parallel, delayed

revio_bam_file_dir = "/hb/groups/cornejo_lab/matt/hla_capture/pacbio/mapped_bam/"
promethion_bam_file_dir = "/hb/groups/cornejo_lab/matt/hla_capture/ont/mapped_bam/"
output_dir = "/hb/groups/cornejo_lab/matt/hla_capture/coverage/"
os.makedirs(output_dir, exist_ok=True)
regions_file = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/mosdepth/mhc_regions.bed"

threads = 4
mapq_threshold = 20

def get_bam_file_paths():
	os.system(f'find {revio_bam_file_dir} -type f -name "*.chr6.bam" | grep -v "bai" > revio_bam_files.txt')
	os.system(f'find {promethion_bam_file_dir} -type f -name "*.chr6.bam" | grep -v "bai" > promethion_bam_files.txt')

	with open("revio_bam_files.txt", "r") as f1, open("promethion_bam_files.txt", "r") as f2:
		bam_file_paths_list = f1.read().splitlines() + f2.read().splitlines()

	os.remove("revio_bam_files.txt")
	os.remove("promethion_bam_files.txt")

	return bam_file_paths_list

def index_bam(bam_file):
	index = "samtools index -b {}".format(bam_file)
	os.system(index)

def run_mosdepth(bam_file):
	filename = os.path.basename(bam_file)
	prefix = os.path.join(output_dir, filename.replace(".bam", ""))
	# --flag 3328 excludes duplicates and secondary/supplementary alignments
	mosdepth = f"mosdepth --flag 3328 --mapq {mapq_threshold} --by {regions_file} --thresholds 20,30 -t {threads} {prefix} {bam_file}"
	os.system(mosdepth)

def main():
	bam_file_paths_list = get_bam_file_paths()
	Parallel(n_jobs=4)(delayed(index_bam)(bam_file) for bam_file in bam_file_paths_list)
	Parallel(n_jobs=6)(delayed(run_mosdepth)(bam_file) for bam_file in bam_file_paths_list)

if __name__ == "__main__":
	main()
