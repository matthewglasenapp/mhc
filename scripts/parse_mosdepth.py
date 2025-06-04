import os
import gzip
import csv
import json
from statistics import mean, stdev

root_dir = "/hb/groups/cornejo_lab/matt/hla_capture/coverage/"
regions_file = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/mosdepth/mhc_regions.bed"
output_dir = "/hb/groups/cornejo_lab/matt/hla_capture/coverage/results/"
biotype_dict_file = "/hb/groups/cornejo_lab/matt/hla_capture/input_data/mosdepth/biotype_dict.json"
os.makedirs(output_dir, exist_ok=True)

threads = 4
# Coordinates for per-base files
start_pos = 28000000
stop_pos = 34000000

per_base_dict = {str(pos): {"revio": [], "promethion": []} for pos in range(start_pos, stop_pos + 1)}

platforms = ["revio", "promethion"]

coverage_dict = {platform: {"gene": {}, "exon": {}} for platform in platforms}

# Ignore HG01891 - failed across all platforms
sample_list = ['HG002', 'HG003', 'HG004', 'HG005', 'HG01106', 'HG01258', 'HG01928', 'HG02055', 'HG02630', 'HG03492', 'HG03579', 'IHW09021', 'IHW09049', 'IHW09071', 'IHW09117', 'IHW09118', 'IHW09122', 'IHW09125', 'IHW09175', 'IHW09198', 'IHW09200', 'IHW09224', 'IHW09245', 'IHW09251', 'IHW09359', 'IHW09364', 'IHW09409', 'NA19240', 'NA20129', 'NA21309', 'NA24694', 'NA24695']

# Original 16 samples
sample_list_hprc = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]

def parse_sample_and_platform(filepath):
	basename = os.path.basename(filepath)

	if "dedup" in basename:
		platform = "revio"
		sample = basename.split(".dedup")[0]
	elif "porechop" in basename:
		platform = "promethion"
		sample = basename.split(".porechop")[0]
	elif "_revio_" in basename:
		platform = "revio"
		sample = basename.split("_revio_")[0]
	elif "_promethion_" in basename:
		platform = "promethion"
		sample = basename.split("_promethion_")[0]
	else:
		raise ValueError(f"Unknown platform in filename: {basename}")

	return sample, platform

def parse_mosdepth_per_base(output_file):
	# 1. Get per-base files
	get_per_base_file_paths = "find {} -type f -name '*per-base.bed.gz*' | grep -v 'csi' > per_base_files".format(root_dir)
	os.system(get_per_base_file_paths)

	with open("per_base_files", "r") as f:
		per_base_file_list = f.read().splitlines()

	os.system("rm per_base_files")

	# 2. Parse per-base files
	for file in per_base_file_list:
		sample, platform = parse_sample_and_platform(file)
		intermediate_output_file = os.path.join(output_dir, f"{sample}_{platform}_per_base.tsv")
		
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
	get_per_base_file_paths = "find {} -type f -name '*_per_base.tsv*' > out_files".format(output_dir)
	os.system(get_per_base_file_paths)

	with open("out_files", "r") as f:
		file_list = f.read().splitlines()

	os.system("rm out_files")
	sorted_file_list = sorted(file_list)
	
	for file in sorted_file_list:
		with open(file, "r") as f:
			sample, platform = parse_sample_and_platform(file)


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
	# Only for the 16 HPRC samples!!
	hprc_indices = [sample_list.index(sample) for sample in sample_list_hprc]
	for platform_name in platforms:
		metric_file = f"{platform_name}_mean_std_depth.tsv"

		with open(output_dir + metric_file, "w") as f1:
			writer = csv.writer(f1, delimiter="\t")

			writer.writerow(["base", "mean_depth", "std_depth"])

			for base, platform_data in per_base_dict.items():
				depths = [platform_data[platform_name][idx] for idx in hprc_indices if idx < len(platform_data[platform_name])]
				avg_depth = mean(depths)
				std_depth = stdev(depths)

				# Write results
				writer.writerow([base, avg_depth, std_depth])

def parse_mosdepth_regions_thresholds(output_json_file):
	with open(biotype_dict_file, "r") as json_file:
		biotype_dict = json.load(json_file)

	for platform in platforms:
		get_regions_file_paths = f"find {root_dir} -type f -name '*.regions.bed.gz' | grep -v 'csi' > regions_files"
		get_thresholds_file_paths = f"find {root_dir} -type f -name '*.thresholds.bed.gz' | grep -v 'csi' > thresholds_files"

		os.system(get_thresholds_file_paths)
		os.system(get_regions_file_paths)

		with open("regions_files", "r") as f1, open("thresholds_files", "r") as f2:
			regions_lines = sorted(f1.read().splitlines())
			thresholds_lines = sorted(f2.read().splitlines())
			print("🧪 Number of regions files found:", len(regions_lines))
			print("🧪 Number of thresholds files found:", len(thresholds_lines))
			file_list = zip(regions_lines, thresholds_lines)

		os.system("rm regions_files")
		os.system("rm thresholds_files")
		
		for regions_file, thresholds_file in file_list:
			print(f"🗂 Processing regions file: {regions_file}")
			print(f"🗂 Processing thresholds file: {thresholds_file}")

			sample_name, platform = parse_sample_and_platform(regions_file)

			if sample_name not in sample_list:
				print(f"⚠️ Skipping unrecognized sample: {sample_name}")
				continue

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
							biotype = biotype_dict[ID]
						elif len(regions_fields[3].split("_")) ==4:
								name = "_".join(regions_fields[3].split("gene_")[1].split("_")[0:2])
								ID = regions_fields[3].split("_")[-1]
								biotype = biotype_dict[ID]
						else:
							ID = regions_fields[3].split("_")[1]
							name = ID
							biotype = biotype_dict[ID]

					elif record_type == "exon":
						name = regions_fields[3]
						ID = name.split("_")[1]
						transcript_ID = name.split("_")[2]
						parent_gene = name.split("_",3)[3]

					if platform not in platforms:
						print(f"⚠️ Platform '{platform}' not recognized")
						continue

					if not ID in coverage_dict[platform][record_type]:
						coverage_dict[platform][record_type][ID] = {
							"name": name,
							"start": start,
							"stop": stop,
							"biotype": biotype if record_type == "gene" else None,
							"transcripts": [] if record_type == "exon" else None,
							"parent_gene": [] if record_type == "exon" else None
						}
					
					coverage_dict[platform][record_type][ID][sample_name] = {
						"coverage_depth": coverage_depth,
						"prop_20x": prop_20x,
						"prop_30x": prop_30x
					}

					print(f"✅ Added {record_type} {ID} for sample={sample_name}, platform={platform}")

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
				output_file = "{}{}_{}_{}.tsv".format(output_dir, platform, record_type, metric)

				with open(output_file, "w") as f:
					writer = csv.writer(f, delimiter="\t")

					if record_type == "exon":
						header = ["ID", "start", "stop", "transcripts", "parent_genes"] + sample_list
					elif record_type == "gene":
						header = ["gene_name", "gene_id", "biotype", "start", "stop"] + sample_list
					writer.writerow(header)

					for record_id, record_info in data.items():
						if record_type == "exon":
							row = [record_id, record_info["start"], record_info["stop"], ",".join(record_info["transcripts"]), ",".join(record_info["parent_gene"])]
						elif record_type == "gene":
							row = [record_info["name"], record_id, record_info["biotype"], record_info["start"], record_info["stop"]]

						row.extend(record_info[sample][metric] for sample in sample_list)
						writer.writerow(row)

def main():
	# parse_mosdepth_per_base(output_dir + "hla_per_base.csv")
	parse_mosdepth_regions_thresholds(output_dir + "coverage_dict.json")
	write_results()

if __name__ == "__main__":
	main()
