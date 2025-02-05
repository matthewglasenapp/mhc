import os
import csv
import subprocess
import multiprocessing

bam_dir = "/hb/groups/cornejo_lab/HLA_hybrid_capture/Pacbio/20240711_Twist-HLA-Panel/HiFiBam/"
output_csv = "barcode_counts.csv"

universal_adapter = "CGAACATGTAGCTGACTCAGGTCAC"
universal_adapter_rc = "GTGACCTGAGTCAGCTACATGTTCG"
me = "AGATGTGTATAAGAGACAG"
me_rc = "CTGTCTCTTATACACATCT"
SMRTbell = "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT"
SMRTbell_rc = "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT"

def process_bam(bam_file):
	get_sample_name = "samtools view -H {} | awk -F'SM:' '/@RG/ {{print $2}}' | cut -f1".format(bam_file)
	sample_name = subprocess.getoutput(get_sample_name).strip()

	total_reads_cmd = "samtools view -c {}".format(bam_file)
	total_reads = int(subprocess.getoutput(total_reads_cmd).strip())

	smrtbell_cmd = "samtools view {} | cut -f10 | grep '{}' | wc -l".format(bam_file, SMRTbell)
	smrtbell_count = int(subprocess.getoutput(smrtbell_cmd).strip())

	smrtbell_rc_cmd = "samtools view {} | cut -f10 | grep '{}' | wc -l".format(bam_file, SMRTbell_rc)
	smrtbell_rc_count = int(subprocess.getoutput(smrtbell_rc_cmd).strip())

	universal_adapter_count_cmd = "samtools view {} | cut -f10 | grep '{}' | wc -l".format(bam_file, universal_adapter)
	universal_adapter_count = int(subprocess.getoutput(universal_adapter_count_cmd).strip())

	universal_adapter_rc_count_cmd = "samtools view {} | cut -f10 | grep '{}' | wc -l".format(bam_file, universal_adapter_rc)
	universal_adapter_rc_count = int(subprocess.getoutput(universal_adapter_rc_count_cmd).strip())

	me_count_cmd = "samtools view {} | cut -f10 | grep '{}' | wc -l".format(bam_file, me)
	me_count = int(subprocess.getoutput(me_count_cmd).strip())

	me_rc_count_cmd = "samtools view {} | cut -f10 | grep '{}' | wc -l".format(bam_file, me_rc)
	me_rc_count = int(subprocess.getoutput(me_rc_count_cmd).strip())

	return [sample_name, total_reads, smrtbell_count, smrtbell_rc_count, universal_adapter_count, universal_adapter_rc_count, me_count, me_rc_count]

def count_barcodes():
	bam_files = [os.path.join(bam_dir, f) for f in os.listdir(bam_dir) if f.endswith(".bam")]

	# Use available CPU cores
	num_processes = min(len(bam_files), multiprocessing.cpu_count())

	with open(output_csv, "w", newline="") as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(["Sample", "TotalReads", "SMRTbell", "SMRTbell_RC", "UniversalAdapter", "UniversalAdapter_RC", "ME", "ME_RC"])

		# Parallel execution
		with multiprocessing.Pool(processes=num_processes) as pool:
			results = pool.map(process_bam, bam_files)

		writer.writerows(results)

def main():
	count_barcodes()

if __name__ == "__main__":
	main()
