import os
import subprocess
import shutil

max_threads = 24

output_directory = "/hb/scratch/mglasena/test_ont/"
barcode_config = "/hb/scratch/mglasena/test_ont/sample_barcode_arrangement.txt"
barcode_file = "/hb/scratch/mglasena/test_ont/sample_barcodes.fa"
basecalled_reads = "/hb/groups/cornejo_lab/HLA_hybrid_capture/06_25_24_R1041_LIG_Cornejo_EXP26/Cornejo/06_25_24_R1041_LIG_Cornejo_EXP26_1_drd0.7.2_sup5.0.0.bam"

demux_prefix = "b2f9e1ada541ad6c3b470699dfbdd70ff26e092f_"

barcode_sample_dict = {
'MY_CUSTOM_KIT_barcode01': 'HG002', 'MY_CUSTOM_KIT_barcode02': 'HG003', 'MY_CUSTOM_KIT_barcode03': 'HG004', 'MY_CUSTOM_KIT_barcode04': 'HG005', 'MY_CUSTOM_KIT_barcode05': 'NA24694', 'MY_CUSTOM_KIT_barcode06': 'NA24695', 'MY_CUSTOM_KIT_barcode07': 'HG01106', 'MY_CUSTOM_KIT_barcode08': 'HG01258', 'MY_CUSTOM_KIT_barcode09': 'HG01891', 'MY_CUSTOM_KIT_barcode10': 'HG01928', 'MY_CUSTOM_KIT_barcode11': 'HG02055', 'MY_CUSTOM_KIT_barcode12': 'HG02630', 'MY_CUSTOM_KIT_barcode13': 'HG03579', 'MY_CUSTOM_KIT_barcode14': 'NA19240', 'MY_CUSTOM_KIT_barcode15': 'NA20129', 'MY_CUSTOM_KIT_barcode16': 'NA21309', 'MY_CUSTOM_KIT_barcode17': 'HG03492', 'MY_CUSTOM_KIT_barcode18': 'IHW09071', 'MY_CUSTOM_KIT_barcode19': 'IHW09021', 'MY_CUSTOM_KIT_barcode20': 'IHW09175', 'MY_CUSTOM_KIT_barcode21': 'IHW09049', 'MY_CUSTOM_KIT_barcode22': 'IHW09117', 'MY_CUSTOM_KIT_barcode23': 'IHW09118', 'MY_CUSTOM_KIT_barcode24': 'IHW09122', 'MY_CUSTOM_KIT_barcode25': 'IHW09125', 'MY_CUSTOM_KIT_barcode26': 'IHW09251', 'MY_CUSTOM_KIT_barcode27': 'IHW09359', 'MY_CUSTOM_KIT_barcode28': 'IHW09364', 'MY_CUSTOM_KIT_barcode29': 'IHW09245', 'MY_CUSTOM_KIT_barcode30': 'IHW09409', 'MY_CUSTOM_KIT_barcode31': 'IHW09198', 'MY_CUSTOM_KIT_barcode32': 'IHW09200', 'MY_CUSTOM_KIT_barcode33': 'IHW09224'
}

def run_dorado():
    dorado_cmd = "dorado demux -o {output_dir} --emit-summary --barcode-arrangement {barcode_config} --barcode-sequences {barcode_sequences} --kit-name MY_CUSTOM_KIT --threads {threads} {basecalled_reads}".format(output_dir = output_directory, barcode_config = barcode_config, barcode_sequences = barcode_file, threads = max_threads, basecalled_reads = basecalled_reads)
   
    subprocess.run(dorado_cmd, shell=True, check=True)

def rename_demux_bams():
    for key, value in barcode_sample_dict.items():
        input_file = os.path.join(output_directory, demux_prefix + key + ".bam")
        output_file = os.path.join(output_directory, value + ".bam")
        shutil.copy(input_file, output_file)

def main():
    run_dorado()
    rename_demux_bams()

if __name__ == "__main__":
    main()
