import os
import subprocess

demux_dir = os.getcwd()
demux_prefix = "b2f9e1ada541ad6c3b470699dfbdd70ff26e092f_"

barcode_sample_dict = {
'MY_CUSTOM_KIT_barcode01': 'HG002', 'MY_CUSTOM_KIT_barcode02': 'HG003', 'MY_CUSTOM_KIT_barcode03': 'HG004', 'MY_CUSTOM_KIT_barcode04': 'HG005', 'MY_CUSTOM_KIT_barcode05': 'NA24694', 'MY_CUSTOM_KIT_barcode06': 'NA24695', 'MY_CUSTOM_KIT_barcode07': 'HG01106', 'MY_CUSTOM_KIT_barcode08': 'HG01258', 'MY_CUSTOM_KIT_barcode09': 'HG01891', 'MY_CUSTOM_KIT_barcode10': 'HG01928', 'MY_CUSTOM_KIT_barcode11': 'HG02055', 'MY_CUSTOM_KIT_barcode12': 'HG02630', 'MY_CUSTOM_KIT_barcode13': 'HG03579', 'MY_CUSTOM_KIT_barcode14': 'NA19240', 'MY_CUSTOM_KIT_barcode15': 'NA20129', 'MY_CUSTOM_KIT_barcode16': 'NA21309', 'MY_CUSTOM_KIT_barcode17': 'HG03492', 'MY_CUSTOM_KIT_barcode18': 'IHW09071', 'MY_CUSTOM_KIT_barcode19': 'IHW09021', 'MY_CUSTOM_KIT_barcode20': 'IHW09175', 'MY_CUSTOM_KIT_barcode21': 'IHW09049', 'MY_CUSTOM_KIT_barcode22': 'IHW09117', 'MY_CUSTOM_KIT_barcode23': 'IHW09118', 'MY_CUSTOM_KIT_barcode24': 'IHW09122', 'MY_CUSTOM_KIT_barcode25': 'IHW09125', 'MY_CUSTOM_KIT_barcode26': 'IHW09251', 'MY_CUSTOM_KIT_barcode27': 'IHW09359', 'MY_CUSTOM_KIT_barcode28': 'IHW09364', 'MY_CUSTOM_KIT_barcode29': 'IHW09245', 'MY_CUSTOM_KIT_barcode30': 'IHW09409', 'MY_CUSTOM_KIT_barcode31': 'IHW09198', 'MY_CUSTOM_KIT_barcode32': 'IHW09200', 'MY_CUSTOM_KIT_barcode33': 'IHW09224'
}

def rename_demux_bams():
    for key, value in barcode_sample_dict.items():
        input_file = os.path.join(demux_dir, demux_prefix + key + ".bam")
        output_file = os.path.join(demux_dir, value + ".bam")
        rename_cmd = "cp {input} {output}".format(input = input_file, output = output_file)
        subprocess.run(rename_cmd, shell=True, check=True)

def main():
    rename_demux_bams()

if __name__ == "__main__":
    main()
