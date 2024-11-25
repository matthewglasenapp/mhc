import json
from statistics import mean

input_json_file = "/Users/matt/Desktop/coverage_dict.json"
targeted_genes_file = "/Users/matt/Documents/GitHub/mhc/bed_files/targeted_genes.bed"

# Ignore HG01891 - failed across all platforms
#sample_list = ['HG002', 'HG003', 'HG004', 'HG005', 'HG01106', 'HG01258', 'HG01928', 'HG02055', 'HG02630', 'HG03492', 'HG03579', 'IHW09021', 'IHW09049', 'IHW09071', 'IHW09117', 'IHW09118', 'IHW09122', 'IHW09125', 'IHW09175', 'IHW09198', 'IHW09200', 'IHW09224', 'IHW09245', 'IHW09251', 'IHW09359', 'IHW09364', 'IHW09409', 'NA19240', 'NA20129', 'NA21309', 'NA24694', 'NA24695']
sample_list = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]

# Original 16 samples
hprc_samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]
ihw_samples = ['IHW09021', 'IHW09049', 'IHW09071', 'IHW09117', 'IHW09118', 'IHW09122', 'IHW09125', 'IHW09175', 'IHW09198', 'IHW09200', 'IHW09224', 'IHW09245', 'IHW09251', 'IHW09359', 'IHW09364', 'IHW09409']

twist_prebuilt_panel = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5']

biotype_dict_file = "/Users/matt/Documents/GitHub/mhc/bed_files/biotype_dict.json"

depth_30x_pacbio = []
prop_20x_pacbio = []
prop_30x_pacbio = []

depth_threshold = 20

class_I = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-J", "HLA-L", "HLA-P", "HLA-S", "HLA-V", "HLA-W"]

class_II = ["HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DRB1", "HLA-DRB5", "HLA-DRB6", "HLA-DRB9"]

class_III = ["ABHD16A", "AGER", "AIF1", "AGPAT1", "APOM", "ATP6V1G2", "ATP6V1G2-DDX39B", "ATF6B", "BAG6", "BTNL2", "C2", "C4A", "C4B", "C6orf47", "CFB", "CLIC1", "CSNK2B", "CYP21A2", "DDAH2", "DDX39B", "DXO", "EHMT2", "EGFL8", "ENSG00000244255", "ENSG00000263020", "ENSG00000285085", "ENSG00000289282", "ENSG00000291302", "FKBPL", "GPSM3", "GPANK1", "HSPA1A", "HSPA1B", "HSPA1L", "LST1", "LSM2", "LTA", "LTB", "LY6G5B", "LY6G5C", "LY6G6C", "LY6G6D", "LY6G6F", "LY6G6F-LY6G6D", "MCCD1", "MPIG6B", "MSH5", "MSH5-SAPCD1", "NCR3", "NEU1", "NELFE", "NFKBIL1", "NOTCH4", "PBX2", "PPT2", "PPT2-EGFL8", "PRRC2A", "PRRT1", "RNF5", "SAPCD1", "SKIC2", "SLC44A4", "STK19", "TNF", "TNXB", "TSBP1", "VARS1", "VWA7", "ZBTB12"]

def process_coverage_dict():
    global biotype_dict
    with open(biotype_dict_file, "r") as json_file1:
        biotype_dict = json.load(json_file1)

    with open(input_json_file, "r") as json_file2:
        coverage_dict = json.load(json_file2)

    global pacbio_gene_dict 
    pacbio_gene_dict = coverage_dict['pacbio']['gene']

    for gene, data in pacbio_gene_dict.items():
        depth_list = []
        prop_20x_list = []
        prop_30x_list = []
        for sample in sample_list:
            depth_list.append(data[sample]['coverage_depth'])
            prop_20x_list.append(data[sample]['prop_20x'])
            prop_30x_list.append(data[sample]['prop_30x'])
        
        if min(depth_list) > depth_threshold:
            depth_30x_pacbio.append(gene)
        
        else:
            gene_name = data['name']
            failed_samples = dict()
            
            for sample in sample_list:
                if data[sample]['coverage_depth'] < depth_threshold:
                    failed_samples[sample] = data[sample]['coverage_depth']
            
            #print("{} ({}) had {} samples with mean coverage depth below 30X".format(gene_name, gene, len(failed_samples)))
            #for key,value in failed_samples.items():
                #print(key,value)

        if min(prop_20x_list) >= 1:
            prop_20x_pacbio.append(gene)
        if min(prop_30x_list) >= 1:
            prop_30x_pacbio.append(gene)

def compare_groups():
    class_I_depth = []
    class_II_depth = []
    class_III_depth = []

    for gene in depth_30x_pacbio:
        depth_values = []
        for sample in sample_list:
            depth_values.append(pacbio_gene_dict[gene][sample]['coverage_depth'])
        mean_depth = mean(depth_values)
        if pacbio_gene_dict[gene]['name'] in class_I:
            class_I_depth.append(mean_depth)
        elif pacbio_gene_dict[gene]['name'] in class_II:
            class_II_depth.append(mean_depth)
        elif pacbio_gene_dict[gene]['name'] in class_III:
            class_III_depth.append(mean_depth)
        else:
            if pacbio_gene_dict[gene]["name"] == "missing":
                print(gene)
            else:
                print(pacbio_gene_dict[gene]["name"])

    c1d = mean(class_I_depth)
    c2d = mean(class_II_depth)
    c3d = mean(class_III_depth)

    print("Mean Class I coverage depth: {:.0f}".format(c1d))
    print("Mean Class II coverage depth: {:.0f}".format(c2d))
    print("Mean Class III coverage depth: {:.0f}".format(c3d))
    print("Enrichment of Class I over class III: {:.1f}X".format(c1d/c3d))
    print("Enrichment of Class II over Class III: {:.1f}X".format(c2d/c3d))
    print("Enrichment of Class I/II over Class III: {:.1f}".format(((c1d+c2d)/2)/c3d))

def print_output():
    print("{} genes had a mean depth of 30X or more across all samples".format(len(depth_30x_pacbio)))
    print("{} genes had all of their bases covered by at least 20 reads".format(len(prop_20x_pacbio)))
    print("{} genes had all of their bases covered by at least 30 reads".format(len(prop_30x_pacbio)))

    targeted_genes = set([gene for item in open(targeted_genes_file, "r").read().splitlines() for gene in item.split("\t")[3].split(";")])
    targeted_genes.update(twist_prebuilt_panel)
    # Ignore HLA-DRB3 and HLA-DRB4 for now because they're not in the hg 38 gtf/gff
    targeted_genes.difference_update({"HLA-DRB3", "HLA-DRB4"})

    # Define captured genes as 30X depth across all samples
    captured_genes = []
    for gene in depth_30x_pacbio:
        if pacbio_gene_dict[gene]["name"] == "missing":
            captured_genes.append(gene)
        else:
            gene_name = pacbio_gene_dict[gene]["name"]
            captured_genes.append(gene_name)

    if targeted_genes.issubset(captured_genes):
        print("All {} targeted genes had a mean coverage depth greater than 20X".format(len(targeted_genes)))
    else:
        overlap = len(set(captured_genes).intersection(targeted_genes))
        failed = set(targeted_genes).difference(captured_genes)
        print("{}/{} targeted genes had a mean coverage depth greater than 30X".format(overlap, len(targeted_genes)))
        print("{} failed genes: {}".format(len(failed),failed))

    captured_protein_coding = 0
    for gene in depth_30x_pacbio:
        if biotype_dict[gene] == "protein_coding":
            captured_protein_coding += 1
        else:
            name = pacbio_gene_dict[gene]["name"]
            #if name == "missing":
                #print(gene)
            #else:
                #print(name)
    #print(captured_protein_coding)

def main():
    process_coverage_dict()
    #print_output()
    compare_groups()

if __name__ == "__main__":
    main()