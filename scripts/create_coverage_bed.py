import sys

HG002_per_base = "/hb/scratch/mglasena/MHC/results/HG002_revio_per_base.tsv"
output_file = "/hb/scratch/mglasena/test2/HG002_revio_30X.bed"

depth_threshold = 30
chromosome = "chr6"

def create_30X_bed(per_base_file, output_file):
	with open(per_base_file, "r") as f1, open(output_file, "w") as f2:
		start = None
		last_pos = None

		for line in f1:
			pos, depth = map(int, line.strip().split("\t"))

			if depth >= depth_threshold:
				if start is None:
					start = pos
				last_pos = pos
			
			else:
				if start is not None:
					f2.write(f"{chromosome}\t{start}\t{pos}\n")
					start = None

		# Handle last region
		if start is not None and last_pos is not None:
			f2.write(f"{chromosome}\t{start}\t{last_pos + 1}\n")

	print("BED file written to {output_file}".format(output_file = output_file))

def main():
	create_30X_bed(HG002_per_base, output_file)

if __name__ == "__main__":
	main()

