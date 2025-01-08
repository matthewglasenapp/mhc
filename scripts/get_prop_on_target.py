import os
import pysam
import pybedtools
import statistics

mapped_bam_dir = "/hb/scratch/mglasena/MHC/mapped_bam/revio/"
output_directory = "/hb/scratch/mglasena/MHC/on_target/"
#bed_file = "/hb/scratch/mglasena/MHC/on_target/mhc.bed"
bed_file = "/hb/scratch/mglasena/MHC/on_target/expected_coverage.bed"

class BAMStatistics:
    def __init__(self, sample):
        self.sample = sample
        self.bam_file = mapped_bam_dir + sample + ".hg38.revio_marked_duplicates.bam"
        self.output_directory = output_directory + sample + "/"
        os.system("mkdir -p {}".format(self.output_directory))

        # Initialize variables
        self.total_reads = 0
        self.total_bases = 0
        self.all_read_lengths = []

        self.on_target_bases = 0
        self.total_target_length = 0
        self.bed_read_lengths = []
        self.overlapping_read_count = 0
        self.average_on_target_depth = 0

    def filter_bam(self):
        # Filter BAM for primary alignments, MAPQ ≥ 20, and non-duplicates
        filtered_bam_file = self.output_directory + self.sample + "_filtered.bam"
        filter_bam = "samtools view -h -q 20 -F 0x4 -F 0x100 -F 0x400 -F 0x800 {} | samtools view -b > {}".format(self.bam_file, filtered_bam_file)
        index_filtered_bam = "samtools index {}".format(filtered_bam_file)
        os.system(filter_bam)
        os.system(index_filtered_bam)

        # Calculate total reads and their lengths for the filtered BAM
        with pysam.AlignmentFile(filtered_bam_file, "rb") as bamfile:
            for read in bamfile:
                self.total_reads += 1
                self.total_bases += read.query_alignment_length
                self.all_read_lengths.append(read.query_length)

        # Calculate total target region length
        self.total_target_length = sum(
            int(feature.end) - int(feature.start) for feature in pybedtools.BedTool(bed_file)
        )

    def intersect_bam_with_targets(self):
        # Intersect filtered BAM with target regions
        filtered_bam_file = self.output_directory + self.sample + "_filtered.bam"
        intersect_output_file = self.output_directory + self.sample + "_on_target.bam"
        intersect = "bedtools intersect -abam {} -b {} > {}".format(filtered_bam_file, bed_file, intersect_output_file)
        index = "samtools index {}".format(intersect_output_file)
        os.system(intersect)
        os.system(index)

    def calculate_statistics(self):
        # Calculate statistics for the intersected BAM
        intersected_bam_file = self.output_directory + self.sample + "_on_target.bam"
        with pysam.AlignmentFile(intersected_bam_file, "rb") as bamfile:
            for read in bamfile:
                self.overlapping_read_count += 1
                self.bed_read_lengths.append(read.query_length)
                self.on_target_bases += read.query_alignment_length

        # Calculate average on-target depth
        if self.total_target_length > 0:
            self.average_on_target_depth = self.on_target_bases / self.total_target_length

    def print_results(self):
        # Proportion calculations
        proportion_on_target_bases = self.on_target_bases / self.total_bases if self.total_bases > 0 else 0
        proportion_reads_overlapping = self.overlapping_read_count / self.total_reads if self.total_reads > 0 else 0

        # Median read lengths
        median_all_reads = statistics.median(self.all_read_lengths) if self.all_read_lengths else 0
        median_bed_reads = statistics.median(self.bed_read_lengths) if self.bed_read_lengths else 0

        # Print results
        print("Total Target Region Length: {} bp".format(self.total_target_length))
        print("Total Reads in BAM: {}".format(self.total_reads))
        print("Reads overlapping target region(s): {}".format(self.overlapping_read_count))
        print("Proportion of Reads Mapping to Target Region: {:.6f}".format(proportion_reads_overlapping))
        print("Median Read Length (All Reads): {}".format(median_all_reads))
        print("Median Read Length (BED-Overlapping Reads): {}".format(median_bed_reads))
        print("Total Aligned Bases: {}".format(self.total_bases))
        print("Total On-Target Bases: {}".format(self.on_target_bases))
        print("Proportion On-Target: {:.6f}".format(proportion_on_target_bases))
        print("Average On-Target Depth: {:.2f}X".format(self.average_on_target_depth))

def main():
    sample_name = "HG002"
    
    sample = BAMStatistics(sample_name)
    sample.filter_bam()
    sample.intersect_bam_with_targets()
    sample.calculate_statistics()
    sample.print_results()

if __name__ == "__main__":
    main()


