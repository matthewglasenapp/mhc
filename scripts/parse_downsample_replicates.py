import os
import csv

root_dir = "/hb/scratch/mglasena/downsample_replicates/"
mhc_classes = ["MHC_Class_I", "MHC_Class_II", "MHC_Class_III"]
proportion_retain = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 
					 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1]
num_replicates = 10

combined_csv = os.path.join(root_dir, "all_concordance_results.csv")

with open(combined_csv, "w", newline="") as outfile:
	writer = csv.writer(outfile)
	writer.writerow(["MHC_Class", "Proportion", "Replicate", "Depth", "Variant", "Metric", "Value"])

	for mhc_class in mhc_classes:
		for proportion in proportion_retain:
			for rep in range(num_replicates):
				summary_path = os.path.join(
					root_dir, str(proportion), f"rep{rep}", f"{mhc_class}_rep{rep}_summary.csv"
				)
				if not os.path.exists(summary_path):
					continue

				with open(summary_path, "r") as infile:
					reader = csv.DictReader(infile)
					for row in reader:
						try:
							writer.writerow([
								mhc_class,
								float(proportion),  # Taken from loop, not the row
								int(rep),
								float(row["Depth"]),
								row["Variant"],
								row["Metric"],
								float(row["Value"])
							])
						except (ValueError, KeyError) as e:
							print(f"Skipping row due to error in {summary_path}: {e}")
