import os
import gzip
import csv
import json
from statistics import mean, stdev

root_dir = "/hb/scratch/mglasena/mhc/"

regions_file = "mhc_regions.bed"

threads = 4

# Coordinates for per-base files
start_pos = 28000000
stop_pos = 34000000

#per_base_dict = {str(pos): {"revio": [], "promethion": [], "minion": []} for pos in range(start_pos, stop_pos + 1)}
per_base_dict = {str(pos): {"pacbio": [], "promethion": []} for pos in range(start_pos, stop_pos + 1)}

#platforms = ["minion", "promethion", "revio"]
platforms = ["pacbio", "promethion"]

coverage_dict = {platform: {"gene": {}, "exon": {}} for platform in platforms}

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
		intermediate_output_file = "{}_per_base.tsv".format(prefix)
		
		with gzip.open(file, "rt") as f1, open(intermediate_output_file, "w") as f2:
			for line in f1:
				if not line.startswith("chr6"):
					continue

				fields = line.strip().split("\t")
				chromosome = fields[0]
				window_start = int(fields[1])
				window_stop = int(fields[2])
				depth = int(fields[3])
				
				# Continue if window end is before gene start 
				if window_stop < start_pos:
					continue
				
				sliding_index = max(window_start, start_pos)
				end_index = min(window_stop, stop_pos)
				while sliding_index <= end_index:
					f2.write(str(sliding_index) + "\t" + str(depth) + "\n")
					sliding_index += 1

	# 3. Collate per-base files
	get_per_base_file_paths = "find {} -type f -name '*_per_base.tsv*' > out_files".format(root_dir)
	os.system(get_per_base_file_paths)

	with open("out_files", "r") as f:
		file_list = f.read().splitlines()

	os.system("rm out_files")
	sorted_file_list = sorted(file_list)
	
	for file in sorted_file_list:
		with open(file, "r") as f:
			sample = file.split("/")[-1].split("per_base.tsv")[0].split("_")[0]
			platform = file.split("/")[-1].split("per_base.tsv")[0].split("_")[1]

			for line in f:
				base, depth = line.strip().split("\t")
				per_base_dict[base][platform].append(int(depth))

	# 4. Write output
	with open(output_file,"w") as outfile:
		writer = csv.writer(outfile, delimiter=",")

		header = ["base", "platform"] + sample_list
		writer.writerow(header)

		for platform_name in platforms:
			for position, platform_data in per_base_dict.items():
				coverage_list = platform_data[platform_name]
				data = [position, platform_name] + coverage_list
				writer.writerow(data)

	# 5. Calculate average and standard deviation for depth by base for each platform 
	for platform_name in platforms:
		metric_file = f"{platform_name}_mean_std_depth.tsv"

		with open(metric_file, "w") as f1:
			writer = csv.writer(f1, delimiter="\t")

			writer.writerow(["base", "mean_depth", "std_depth"])

			for base, platform_data in per_base_dict.items():
				depths = platform_data[platform_name]
				avg_depth = mean(depths)
				std_depth = stdev(depths)

				# Write results
				writer.writerow([base, avg_depth, std_depth])

def parse_mosdepth_regions_thresholds(output_json_file):
	for platform in platforms:

		get_regions_file_paths = "find {} -type f -name '*.regions.bed.gz*' | grep {} | grep -v 'csi' > regions_files".format(root_dir, platform)
		get_thresholds_file_paths = "find {} -type f -name '*.thresholds.bed.gz*' | grep {} | grep -v 'csi' > thresholds_files".format(root_dir, platform)
		
		os.system(get_thresholds_file_paths)
		os.system(get_regions_file_paths)

		with open("regions_files", "r") as f1, open("thresholds_files","r") as f2:
			file_list = zip(sorted(f1.read().splitlines()),sorted(f2.read().splitlines()))

		os.system("rm regions_files")
		os.system("rm thresholds_files")
		
		for regions_file, thresholds_file in file_list:
			sample_name = regions_file.split("/")[-1].split(".")[0].split("_")[0]

			with gzip.open(regions_file, "rt") as f1, gzip.open(thresholds_file,"rt") as f2:
				regions = f1.read().splitlines()
				thresholds = f2.read().splitlines()[1:]

				for regions_line, thresholds_line in zip(regions,thresholds):
					regions_fields = regions_line.split("\t")
					record_type = regions_fields[3].split("_")[0]
					start = regions_fields[1]
					stop = regions_fields[2]
					length = int(stop) - int(start)
					coverage_depth = float(regions_fields[4])
					threshold_fields = thresholds_line.split("\t")
					num_20x = int(threshold_fields[4])
					num_30x = int(threshold_fields[5])
					prop_20x = num_20x / length
					prop_30x = num_30x / length
					
					if record_type == "gene":
						if len(regions_fields[3].split("_")) == 3:
							name = regions_fields[3].split("_")[1]
							ID = regions_fields[3].split("_")[2]
						else:
							name = "missing"
							ID = regions_fields[3].split("_")[1]
					
					elif record_type == "exon":
						name = regions_fields[3]
						ID = name.split("_")[1]
						transcript_ID = name.split("_")[2]
						parent_gene = name.split("_",3)[3]

					if not ID in coverage_dict[platform][record_type]:
						coverage_dict[platform][record_type][ID] = {
							"name": name,
							"start": start,
							"stop": stop,
							"transcripts": [] if record_type == "exon" else None,
							"parent_gene": [] if record_type == "exon" else None
						}
					
					coverage_dict[platform][record_type][ID][sample_name] = {
						"coverage_depth": coverage_depth,
						"prop_20x": prop_20x,
						"prop_30x": prop_30x
					}

					current_record = coverage_dict[platform][record_type][ID]
					if record_type == "exon":
						if not transcript_ID in current_record["transcripts"]:
							current_record["transcripts"].append(transcript_ID)
						if not parent_gene in current_record["parent_gene"]:
							current_record["parent_gene"].append(parent_gene)

	with open(output_json_file, "w") as json_file:
		json.dump(coverage_dict, json_file, indent = 4)

def write_results():
	metrics = ["coverage_depth", "prop_20x", "prop_30x"]

	for platform, record_types in coverage_dict.items():
		for record_type, data in record_types.items():
			for metric in metrics:
				output_file = f"{platform}_{record_type}_{metric}.tsv"

				with open(output_file, "w") as f:
					writer = csv.writer(f, delimiter="\t")

					if record_type == "exon":
						header = ["ID", "start", "stop", "transcripts", "parent_genes"] + sample_list
					elif record_type == "gene":
						header = ["gene_name", "gene_id", "start", "stop"] + sample_list
					writer.writerow(header)

					for record_id, record_info in data.items():
						if record_type == "exon":
							row = [record_id, record_info["start"], record_info["stop"], ",".join(record_info["transcripts"]), ",".join(record_info["parent_gene"])]
						elif record_type == "gene":
							row = [record_info["name"], record_id, record_info["start"], record_info["stop"]]

						row.extend(record_info[sample][metric] for sample in sample_list)
						writer.writerow(row)

def main():
	#parse_mosdepth_per_base("hla_per_base.csv")
	parse_mosdepth_regions_thresholds("coverage_dict.json")
	write_results()

if __name__ == "__main__":
	main()
