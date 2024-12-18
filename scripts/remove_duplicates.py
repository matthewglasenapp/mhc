import os
from joblib import Parallel, delayed

bam_file_dir = "/hb/scratch/mglasena/MHC/mapped_bam/"

def get_file_paths_list():
    get_bam_paths_file = 'find {} -type f -name *.bam* | grep -v "bai" > bam_files.txt'.format(bam_file_dir)
    os.system(get_bam_paths_file)

    with open("bam_files.txt", "r") as f:
        bam_file_paths_list = f.read().splitlines()

    os.system("rm bam_files.txt")

    return bam_file_paths_list

def remove_duplicates(bam_file):
    output_file = bam_file.split("_marked_duplicates.bam")[0] + "_rmdup.bam"

    rmdup = "samtools view -h -F 1024 -o {} {}".format(output_file, bam_file)
    os.system(rmdup)
    
def main():
    bam_file_paths_list = get_file_paths_list()
    Parallel(n_jobs=len(bam_file_paths_list))(delayed(remove_duplicates)(bam_file) for bam_file in bam_file_paths_list)

if __name__ == "__main__":
    main()

