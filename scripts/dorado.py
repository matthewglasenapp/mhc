import os
import subprocess
import shutil

max_threads = 24

# set input directory to the current working directory where the script should be run
input_dir = os.getcwd()

# Set the output directory to processed_data/
output_dir = os.path.join(input_dir, "processed_data")
os.makedirs(output_dir, exist_ok=True)

# Use reference fasta with no alternate contigs.
reference_fasta = os.path.join(input_dir, "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")

# Set mapped chr6 reads threshold at which variant calling should not proceed
min_reads_sample = 100

output_directory = "/hb/scratch/mglasena/test_ont/"
barcode_config = "/hb/scratch/mglasena/test_ont/sample_barcode_arrangement.txt"
barcode_file = "/hb/scratch/mglasena/test_ont/sample_barcodes.fa"
basecalled_reads = "/hb/groups/cornejo_lab/HLA_hybrid_capture/06_25_24_R1041_LIG_Cornejo_EXP26/Cornejo/06_25_24_R1041_LIG_Cornejo_EXP26_1_drd0.7.2_sup5.0.0.bam"

demux_prefix = "b2f9e1ada541ad6c3b470699dfbdd70ff26e092f_"

barcode_sample_dict = {
'MY_CUSTOM_KIT_barcode01': 'HG002', 'MY_CUSTOM_KIT_barcode02': 'HG003', 'MY_CUSTOM_KIT_barcode03': 'HG004', 'MY_CUSTOM_KIT_barcode04': 'HG005', 'MY_CUSTOM_KIT_barcode05': 'NA24694', 'MY_CUSTOM_KIT_barcode06': 'NA24695', 'MY_CUSTOM_KIT_barcode07': 'HG01106', 'MY_CUSTOM_KIT_barcode08': 'HG01258', 'MY_CUSTOM_KIT_barcode09': 'HG01891', 'MY_CUSTOM_KIT_barcode10': 'HG01928', 'MY_CUSTOM_KIT_barcode11': 'HG02055', 'MY_CUSTOM_KIT_barcode12': 'HG02630', 'MY_CUSTOM_KIT_barcode13': 'HG03579', 'MY_CUSTOM_KIT_barcode14': 'NA19240', 'MY_CUSTOM_KIT_barcode15': 'NA20129', 'MY_CUSTOM_KIT_barcode16': 'NA21309', 'MY_CUSTOM_KIT_barcode17': 'HG03492', 'MY_CUSTOM_KIT_barcode18': 'IHW09071', 'MY_CUSTOM_KIT_barcode19': 'IHW09021', 'MY_CUSTOM_KIT_barcode20': 'IHW09175', 'MY_CUSTOM_KIT_barcode21': 'IHW09049', 'MY_CUSTOM_KIT_barcode22': 'IHW09117', 'MY_CUSTOM_KIT_barcode23': 'IHW09118', 'MY_CUSTOM_KIT_barcode24': 'IHW09122', 'MY_CUSTOM_KIT_barcode25': 'IHW09125', 'MY_CUSTOM_KIT_barcode26': 'IHW09251', 'MY_CUSTOM_KIT_barcode27': 'IHW09359', 'MY_CUSTOM_KIT_barcode28': 'IHW09364', 'MY_CUSTOM_KIT_barcode29': 'IHW09245', 'MY_CUSTOM_KIT_barcode30': 'IHW09409', 'MY_CUSTOM_KIT_barcode31': 'IHW09198', 'MY_CUSTOM_KIT_barcode32': 'IHW09200', 'MY_CUSTOM_KIT_barcode33': 'IHW09224'
}

# Dictionary of samples to process.
# Sample ID: Read Group String
sample_dict = {
    "HG002" : "@RG\tID:m84039_240622_113450_s1\tSM:HG002",
    "HG003" : "@RG\tID:m84039_240622_113450_s1\tSM:HG003",
    "HG004" : "@RG\tID:m84039_240622_113450_s1\tSM:HG004",
    "HG005" : "@RG\tID:m84039_240622_113450_s1\tSM:HG005",
    "HG01106" : "@RG\tID:m84039_240622_113450_s1\tSM:HG01106",
    "HG01258" : "@RG\tID:m84039_240622_113450_s1\tSM:HG01258",
    "HG01891" : "@RG\tID:m84039_240622_113450_s1\tSM:HG01891",
    "HG01928" : "@RG\tID:m84039_240622_113450_s1\tSM:HG01928",
    "HG02055" : "@RG\tID:m84039_240622_113450_s1\tSM:HG02055",
    "HG02630" : "@RG\tID:m84039_240622_113450_s1\tSM:HG02630",
    "HG03492" : "@RG\tID:m84039_240622_113450_s1\tSM:HG03492",
    "HG03579" : "@RG\tID:m84039_240622_113450_s1\tSM:HG03579",
    "IHW09021" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09021",
    "IHW09049" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09049",
    "IHW09071" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09071",
    "IHW09117" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09117",
    "IHW09118" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09118",
    "IHW09122" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09122",
    "IHW09125" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09125",
    "IHW09175" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09175",
    "IHW09198" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09198",
    "IHW09200" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09200",
    "IHW09224" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09224",
    "IHW09245" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09245",
    "IHW09251" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09251",
    "IHW09359" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09359",
    "IHW09364" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09364",
    "IHW09409" : "@RG\tID:m84039_240622_113450_s1\tSM:IHW09409",
    "NA19240" : "@RG\tID:m84039_240622_113450_s1\tSM:NA19240",
    "NA20129" : "@RG\tID:m84039_240622_113450_s1\tSM:NA20129",
    "NA21309" : "@RG\tID:m84039_240622_113450_s1\tSM:NA21309",
    "NA24694" : "@RG\tID:m84039_240622_113450_s1\tSM:NA24694",
    "NA24695" : "@RG\tID:m84039_240622_113450_s1\tSM:NA24695"
}

# Ensure all required tools are installed and executable
def check_required_commands():    
    print("Checking the installation status of the required bioinformatics tools!")

    required_commands = [
        "bam2fastq",
        "bcftools",
        "bgzip",
        "dorado",
        "fastqc",
        "pigz",
        "porechop_abi",
        "samtools",
        "singularity",
        "tabix",
        "whatshap"
    ]

    missing_commands = []
    for command in required_commands:
        if shutil.which(command) is None:
            missing_commands.append(command)
    if len(missing_commands) != 0:
        print("Error: Missing the following commands: {}".format(", ".join(missing_commands)))
        sys.exit(1)
    else:
        print("All tools required are installed!")
        print("\n\n")

def run_dorado():
    dorado_cmd = "dorado demux -o {output_dir} --emit-summary --barcode-arrangement {barcode_config} --barcode-sequences {barcode_sequences} --kit-name MY_CUSTOM_KIT --threads {threads} {basecalled_reads}".format(output_dir = output_directory, barcode_config = barcode_config, barcode_sequences = barcode_file, threads = max_threads, basecalled_reads = basecalled_reads)
   
    subprocess.run(dorado_cmd, shell=True, check=True)

def rename_demux_bams():
    for key, value in barcode_sample_dict.items():
        input_file = os.path.join(output_directory, demux_prefix + key + ".bam")
        output_file = os.path.join(output_directory, "raw_bam", value + ".bam")
        shutil.copy(input_file, output_file)

class Samples:
    raw_bam_dir = os.path.join(output_dir, "raw_bam")
    mapped_bam_dir = os.path.join(output_dir, "mapped_bam")
    clair3_dir = os.path.join(output_dir, "clair3_vcf")
    whatshap_phased_vcf_dir = os.path.join(output_dir, "phased_vcf_whatshap")

    def __init__(self, sample_ID, read_group_string):
        self.sample_ID = sample_ID
        self.unmapped_bam = os.path.join(raw_bam_dir, self.sample_ID + ".bam")
        self.read_group_string = read_group_string

        for directory in [Samples.raw_bam_dir, Samples.mapped_bam_dir, Samples.clair3_dir, Samples.whatshap_phased_vcf_dir]:
            os.makedirs(directory, exist_ok=True)
        
        print(f"Processing Sample {sample_ID}!")
        print("\n\n")
    
    # Convert BAM file of unmapped HiFi (ccs) reads to FASTQ format for marking duplicates and trimming adapters
    def convert_bam_to_fastq(self):
        print("Converting raw reads to fastq format using pbtk bam2fastq!")
        print("bam2fastq input file: {}".format(self.unmapped_bam))
        
        os.chdir(Samples.fastq_raw_dir)
        
        bam2fastq_cmd = "bam2fastq -j {threads} {input_file} -o {output_prefix}".format(threads = max_threads, input_file = self.unmapped_bam, output_prefix = self.sample_ID)
        
        subprocess.run(bam2fastq_cmd, shell=True, check=True)
                
        print("Raw fastq reads written to: {}".format(os.path.join(Samples.fastq_raw_dir, self.sample_ID + ".fastq.gz")))
        print("\n\n")

    def run_porechop_abi(self):
        input_fastq = ""
        output_fastq = ""
        porechop_cmd = "porechop_abi --ab_initio -i {input_file} -t {threads} -o {output_file} --format fastq.gz".format(input_file = input_fastq, threads = max_threads, output_file = output_fastq)

    def align_to_reference(self):
        pass

    def mark_duplicates(self):
        pass

    def filter_reads(self):
        pass

    def call_variants(self):
        pass

    def call_structural_variants(self):
        pass

    def phase_genotypes_whatshap(self):
        pass


def main():
    # Check that all required tools are installed
    check_required_commands()

    run_dorado()
    rename_demux_bams()

    for sample_ID, sample_read_group_string in sample_dict.items():
        start_time = time.time()
        sample = Samples(sample_ID, sample_read_group_string)
        sample.convert_bam_to_fastq()
        sample.run_porechop_abit()
        sample.align_to_reference()
        sample.mark_duplicates()

        chr6_reads = sample.filter_reads()

        if chr6_reads > min_reads_sample:
            sample.call_variants()
            sample.call_structural_variants()
            sample.phase_genotypes_whatshap()

        else:
            print("Insufficient reads for variant calling")
            print("Sample {sample_id} had {num_reads} reads!".format(sample_id = sample_ID, num_reads = chr6_reads))

if __name__ == "__main__":
    main()
