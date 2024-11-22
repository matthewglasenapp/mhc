import os
import gzip
import csv
from joblib import Parallel, delayed

root_dir = "/hb/scratch/mglasena/hla/mhc/"

regions_file = "mhc_regions.bed"

threads = 4

# Coordinates for per-base files
start_pos = 28000000
stop_pos = 34000000

per_base_dict = {str(pos): {"revio": [], "promethion": [], "minion": []} for pos in range(start_pos, stop_pos + 1)}
coverage_dict_genes = dict()
coverage_dict_exons = dict()

platforms = ["minion", "promethion", "revio"]

sample_list = ['HG002', 'HG003', 'HG004', 'HG005', 'HG01106', 'HG01258', 'HG01891', 'HG01928', 'HG02055', 'HG02630', 'HG03492', 'HG03579', 'IHW09021', 'IHW09049', 'IHW09071', 'IHW09117', 'IHW09118', 'IHW09122', 'IHW09125', 'IHW09175', 'IHW09198', 'IHW09200', 'IHW09224', 'IHW09245', 'IHW09251', 'IHW09359', 'IHW09364', 'IHW09409', 'NA19240', 'NA20129', 'NA21309', 'NA24694', 'NA24695']

def parse_mosdepth_per_base(output_file):
	# 1. Get per-base files
	get_per_base_file_paths = "find {} -type f -name '*per-base.bed.gz*' | grep -v 'csi' > per_base_files".format(root_dir)
	os.system(get_per_base_file_paths)

	with open("per_base_files", "r") as f:
		per_base_file_list = f.read().splitlines()

	os.system("rm per_base_files")

	# 2. Parse per-base files
	for file in per_base_file_list:
		prefix = file.split(".per-base.bed.gz")[0].split("/")[-1]
		output_file = "{}_per_base.tsv".format(prefix)
		lines = gzip.open(file, "rt").read().splitlines()
	
		with open(output_file, "w") as f:
			for line in lines:
				chromosome = line.split("\t")[0]

				if chromosome == "chr6":
					window_start = int(line.split("\t")[1])
					window_stop = int(line.split("\t")[2])
					depth = int(line.split("\t")[3])
					
					# Continue if window end is before gene start 
					if window_stop < start_pos:
						continue
					
					# Window begins before gene start and ends after gene start 
					elif window_stop >= start_pos and window_start < start_pos:
						sliding_index = start_pos
						while sliding_index < window_stop:
							f.write(str(sliding_index) + "\t" + str(depth) + "\n")
							sliding_index += 1

					# Window contained entirely within gene start/stop coordinates
					elif window_start >= start_pos and window_stop <= stop_pos:
						sliding_index = window_start
						while sliding_index < window_stop:
							f.write(str(sliding_index) + "\t" + str(depth) + "\n")
							sliding_index += 1

					# Window begins before gene start but ends after gene start 
					elif window_start >= start_pos and window_start <= stop_pos and window_stop >= stop_pos:
						sliding_index = window_start
						while sliding_index <= stop_pos:
							f.write(str(sliding_index) + "\t" + str(depth) + "\n")
							sliding_index += 1

	# 3. Collate per-base files
	get_per_base_file_paths = "find {} -type f -name '*_per_base.tsv*' > out_files".format(root_dir)
	os.system(get_per_base_file_paths)

	with open("out_files", "r") as f:
		file_list = f.read().splitlines()

	os.system("rm out_files")
	sorted_file_list = sorted(file_list)
	
	for file in sorted_file_list:
		lines = open(file,"r").read().splitlines()
		sample = file.split("/")[-1].split("per_base.tsv")[0].split("_")[0]
		platform = file.split("/")[-1].split("per_base.tsv")[0].split("_")[1]

		for line in lines:
			base = str(line.split("\t")[0])
			depth = str(line.split("\t")[1])
			per_base_dict[base][platform].append(depth)

	# 4. Write output
	outfile = open(output_file,"w")
	writer = csv.writer(out_file, delimiter=",")

	header = ["base", "platform"] + sample_list
	writer.writerow(header)

	platforms = ['revio', 'promethion', 'minion']

	for platform_name in platforms:
		for position, platform_data in per_base_dict.items():
			coverage_list = platform_data[platform_name]
			data = [position, platform_name] + coverage_list
			writer.writerow(data)

	out_file.close()

# Need to account for platform!
def parse_mosdepth_regions_thresholds():
	for platform in platforms:
		coverage_dict_genes[platform] = dict()
		coverage_dict_exons[platform] = dict()

		get_regions_file_paths = "find {} -type f -name '*.regions.bed.gz*' | grep {} | grep -v 'csi' > regions_files".format(root_dir, platform)
		os.system(get_regions_file_paths)

		get_thresholds_file_paths = "find {} -type f -name '*.thresholds.bed.gz*' | grep {} | grep -v 'csi' > thresholds_files".format(root_dir, platform)
		os.system(get_thresholds_file_paths)

		with open("regions_files", "r") as f1, open("thresholds_files","r") as f2:
			file_list = zip(sorted(f1.read().splitlines()),sorted(f2.read().splitlines()))

		os.system("rm regions_files")
		os.system("rm thresholds_files")
		
		for file_pair in list(file_list):
			regions_file = file_pair[0]
			thresholds_file = file_pair[1]

			sample_name = regions_file.split("/")[-1].split(".")[0]

			with gzip.open(regions_file, "rt") as f1, gzip.open(thresholds_file,"rt") as f2:
				records = f1.read().splitlines()
				thresholds = f2.read().splitlines()[1:]

				zipped_list = list(zip(records,thresholds))

				for record in zipped_list:
					# If record is gene
					record_type = record[0].split("\t")[3].split("_")[0]
					
					if record_type == "gene":
						if len(record[0].split("\t")[3].split("_")) == 3:
							name = record[0].split("\t")[3].split("_")[1]
							id = record[0].split("\t")[3].split("_")[2]
						else:
							name = "missing"
							id = record[0].split("\t")[3].split("_")[1]
					
					elif record_type == "exon":
						name = record[0].split("\t")[3]
						id = name.split("_")[1]
						transcript_id = name.split("_")[2]

					start = record[0].split("\t")[1]
					stop = record[0].split("\t")[2]
					length = int(stop) - int(start)
					coverage_depth = float(record[0].split("\t")[4])
					num_20x = int(record[1].split("\t")[4])
					num_30x = int(record[1].split("\t")[5])
					prop_20x = num_20x / length
					prop_30x = num_30x / length
					
					coverage_dict_genes[platform][id] = dict()
					coverage_dict_exons[platform][id] = dict()

					if record_type == "gene":
						coverage_dict_genes[platform][id][sample_name] = [gene_name, start, stop, coverage_depth, prop_20x, prop_30x]

					elif record_type == "exon":
						if not exon_id in coverage_depth_dict[platform][id][sample_name]:
							coverage_dict_exons[platform][id][sample_name] = [exon_name, [transcript_id], [gene_id], start, stop, coverage_depth, prop_20x, prop_30x]

						else:
							coverage_depth_dict[platform][id][sample_name][1].append(transcript_id)
							coverage_depth_dict[platform][id][sample_name][2].append(id)

def main():
	#parse_mosdepth_per_base("hla_per_base.csv")
	#parse_mosdepth_regions_thresholds()

	#Write CSV and JSon

if __name__ == "__main__":
	main()