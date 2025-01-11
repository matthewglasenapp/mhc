import os
import pysam
import pybedtools
import statistics
import csv

platforms = ["revio", "promethion"]

mapped_bam_dir = "/hb/scratch/mglasena/MHC/mapped_bam/"
output_directory = "/hb/scratch/mglasena/MHC/on_target/"
bed_file = "/hb/scratch/mglasena/MHC/on_target/expected_coverage.bed"

sample_list = ['HG002', 'HG003', 'HG004', 'HG005', 'HG01106', 'HG01258', 'HG01928', 'HG02055', 'HG02630', 'HG03492', 'HG03579', 'NA19240', 'NA20129', 'NA21309', 'NA24694', 'NA24695', 'HG01891', 'IHW09021', 'IHW09049', 'IHW09071', 'IHW09117', 'IHW09118', 'IHW09122', 'IHW09125', 'IHW09175', 'IHW09198', 'IHW09200', 'IHW09224', 'IHW09245', 'IHW09251', 'IHW09359', 'IHW09364', 'IHW09409']

data_dict = {sample: {} for sample in sample_list}

def calculate_total_target_length(bed_file):
    total_target_length = 0
    for feature in pybedtools.BedTool(bed_file):
        feature_length = int(feature.end) - int(feature.start)
        total_target_length += feature_length
    return total_target_length

class BAMStatistics:
    total_target_length = 0

    def __init__(self, sample, platform):
        self.sample = sample
        self.platform = platform
        self.bam_file = mapped_bam_dir + platform + "/" + sample + ".hg38." + platform + "_marked_duplicates.bam"
        self.output_directory = output_directory + platform + "/" + sample + "/"
        os.system("mkdir -p {}".format(self.output_directory))

        # Initialize variables
        self.total_reads = 0
        self.total_bases = 0
        self.all_read_lengths = []
        self.on_target_bases = 0
        self.bed_read_lengths = []
        self.overlapping_read_count = 0
        self.average_on_target_depth = 0
        self.reference_coverage_bases = 0

    def filter_bam(self):
        filtered_bam_file = self.output_directory + self.sample + "_filtered.bam"
        filter_bam = "samtools view -h -q 20 -F 0x4 -F 0x100 -F 0x400 -F 0x800 {} | samtools view -b > {}".format(self.bam_file, filtered_bam_file)
        os.system(filter_bam)
        os.system("samtools index {}".format(filtered_bam_file))

        # Process the QC-filtered BAM file
        with pysam.AlignmentFile(filtered_bam_file, "rb") as bamfile:
            for read in bamfile:
                self.total_reads += 1
                aligned_read_bases = 0

                for operation, length in read.cigartuples:
                    # M (match)
                    if operation == 0:
                        aligned_read_bases += length

                    # I (insertion)
                    elif operation == 1:
                        aligned_read_bases += length

                # Update total bases for on-target calculations
                self.total_bases += aligned_read_bases

                self.all_read_lengths.append(read.query_length)



    def intersect_bam_with_targets(self):
        filtered_bam_file = self.output_directory + self.sample + "_filtered.bam"
        intersect_output_file = self.output_directory + self.sample + "_on_target.bam"

        intersect = "bedtools intersect -abam {} -b {} > {}".format(filtered_bam_file, bed_file, intersect_output_file)
        index = "samtools index {}".format(intersect_output_file)
        os.system(intersect)
        os.system(index)

    def calculate_statistics(self):
        intersected_bam_file = self.output_directory + self.sample + "_on_target.bam"

        with pysam.AlignmentFile(intersected_bam_file, "rb") as bamfile:
            for read in bamfile:
                self.overlapping_read_count += 1
                aligned_bases = 0
                coverage_bases = 0 

                for operation, length in read.cigartuples:
                    # M (match)
                    if operation == 0:
                        aligned_bases += length
                        coverage_bases += length

                    # I (insertion)
                    elif operation == 1:
                        aligned_bases += length

                    # D (deletion)
                    elif operation ==2:
                        coverage_bases += length

                self.bed_read_lengths.append(read.query_length)
                self.on_target_bases += aligned_bases
                self.reference_coverage_bases += coverage_bases

        # Calculate the average on-target depth
        self.average_on_target_depth = self.reference_coverage_bases / BAMStatistics.total_target_length

        proportion_on_target_bases = self.on_target_bases / self.total_bases
        proportion_reads_overlapping = self.overlapping_read_count / self.total_reads
        median_all_reads = statistics.median(self.all_read_lengths)
        mean_all_reads = statistics.mean(self.all_read_lengths)
        median_bed_reads = statistics.median(self.bed_read_lengths)
        mean_bed_reads = statistics.mean(self.bed_read_lengths)

        # Populate the results dictionary for this sample
        data_dict[self.sample] = {
            "Total Reads": self.total_reads,
            "Reads Overlapping Target Region": self.overlapping_read_count,
            "Proportion of Reads Mapping to Target Region": proportion_reads_overlapping,
            "Median Read Length (All Reads)": median_all_reads,
            "Mean Read Length (All Reads)": mean_all_reads,
            "Median Read Length (BED-Overlapping Reads)": median_bed_reads,
            "Mean Read Length (BED-Overlapping Reads)": mean_bed_reads,
            "Total Aligned Bases": self.total_bases,
            "On-Target Bases": self.on_target_bases,
            "Proportion On-Target": proportion_on_target_bases,
            "Average On-Target Depth": self.average_on_target_depth
        }

def write_results(csv_file):
    # Write the results to a CSV file
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([
            "Sample", "Total Reads", "Reads Overlapping Target Region",
            "Proportion of Reads Mapping to Target Region", "Median Read Length (All Reads)",
            "Mean Read Length (All Reads)", "Median Read Length (BED-Overlapping Reads)",
            "Mean Read Length (BED-Overlapping Reads)", "Total Aligned Bases", "On-Target Bases",
            "Proportion On-Target", "Average On-Target Depth"
        ])
        for sample, stats in data_dict.items():
            writer.writerow([
                sample,
                stats["Total Reads"],
                stats["Reads Overlapping Target Region"],
                stats["Proportion of Reads Mapping to Target Region"],
                stats["Median Read Length (All Reads)"],
                stats["Mean Read Length (All Reads)"],
                stats["Median Read Length (BED-Overlapping Reads)"],
                stats["Mean Read Length (BED-Overlapping Reads)"],
                stats["Total Aligned Bases"],
                stats["On-Target Bases"],
                stats["Proportion On-Target"],
                stats["Average On-Target Depth"]
            ])

def main():
    BAMStatistics.total_target_length = calculate_total_target_length(bed_file)

    for platform in platforms:
        csv_file = output_directory + platform + "/results.csv"

        for sample_name in sample_list:
            sample = BAMStatistics(sample_name, platform)
            sample.filter_bam()
            sample.intersect_bam_with_targets()
            sample.calculate_statistics()

        write_results(csv_file)

if __name__ == "__main__":
    main()
