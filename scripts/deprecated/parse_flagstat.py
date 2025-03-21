import os
import csv

flagstat_dir = "/hb/scratch/mglasena/MHC/mapped_bam/"
output_dir = "/hb/scratch/mglasena/MHC/results/"

samples = ['HG002', 'HG003', 'HG004', 'HG005', 'HG01106', 'HG01258', 'HG01891', 'HG01928', 'HG02055', 'HG02630', 'HG03492', 'HG03579', 'IHW09021', 'IHW09049', 'IHW09071', 'IHW09117', 'IHW09118', 'IHW09122', 'IHW09125', 'IHW09175', 'IHW09198', 'IHW09200', 'IHW09224', 'IHW09245', 'IHW09251', 'IHW09359', 'IHW09364', 'IHW09409', 'NA19240', 'NA20129', 'NA21309', 'NA24694', 'NA24695']

flagstat_output_dict = {sample: {'promethion': [], 'revio': []} for sample in samples}

def create_flagstat_output_dict():
	get_flagstat_files = "find {} -type f -name *.tsv* > flagstat_files.txt".format(flagstat_dir)
	os.system(get_flagstat_files)

	with open("flagstat_files.txt","r")as f:
		flagstat_files = f.read().splitlines()

	os.system("rm flagstat_files.txt")

	for file in flagstat_files:
		sample = file.split("/")[-1].split("_")[0]
		platform = file.split("/")[-2]
		with open(file,"r") as f:
			lines = f.read().splitlines()
			total = int(lines[0].split("\t")[0])
			duplicates = int(lines[4].split("\t")[0])
			percent_duplicates = (duplicates / float(total)) * 100
			primary_mapped = int(lines[8].split("\t")[0])
			primary_mapped_percent = float(lines[9].split("%")[0])
			flagstat_output_dict[sample][platform] = [total, primary_mapped, primary_mapped_percent, duplicates, percent_duplicates]

def write_results():
    # Define output filenames
    promethion_file = "promethion_flagstat_results.csv"
    revio_file = "revio_flagstat_results.csv"

    # Define column headers
    headers = ["sample", "total", "primary_mapped", "primary_mapped_percent", "duplicates", "percent_duplicates"]

    # Write Promethion results
    with open(output_dir + promethion_file, "w", newline="") as prom_file:
        writer = csv.writer(prom_file)
        writer.writerow(headers)  # Write header row

        for sample in samples:
            if flagstat_output_dict[sample]['promethion']:
                writer.writerow([sample] + flagstat_output_dict[sample]['promethion'])

    # Write Revio results
    with open(output_dir + revio_file, "w", newline="") as rev_file:
        writer = csv.writer(rev_file)
        writer.writerow(headers)  # Write header row

        for sample in samples:
            if flagstat_output_dict[sample]['revio']:
                writer.writerow([sample] + flagstat_output_dict[sample]['revio'])

def main():
	create_flagstat_output_dict()
	write_results()

if __name__ == "__main__":
	main()