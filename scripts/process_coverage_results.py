import json
from statistics import mean

input_json_file = "/Users/matt/Desktop/coverage_dict.json"
targeted_genes_file = "/Users/matt/Documents/GitHub/mhc/bed_files/targeted_genes.bed"
biotype_dict_file = "/Users/matt/Documents/GitHub/mhc/bed_files/biotype_dict.json"

# Ignore HG01891 - failed across all platforms
# Ignore IHW09118, IHW09245, IHW09364, IHW09071  - low coverage
#sample_list = ['HG002', 'HG003', 'HG004', 'HG005', 'HG01106', 'HG01258', 'HG01928', 'HG02055', 'HG02630', 'HG03492', 'HG03579', 'IHW09021', 'IHW09049', 'IHW09117', 'IHW09122', 'IHW09125', 'IHW09175', 'IHW09198', 'IHW09200', 'IHW09224', 'IHW09251', 'IHW09359', 'IHW09409', 'NA19240', 'NA20129', 'NA21309', 'NA24694', 'NA24695']

# Full sample _list
sample_list = ['HG002', 'HG003', 'HG004', 'HG005', 'HG01106', 'HG01258', 'HG01928', 'HG02055', 'HG02630', 'HG03492', 'HG03579', 'NA19240', 'NA20129', 'NA21309', 'NA24694', 'NA24695', 'IHW09021', 'IHW09049', 'IHW09071', 'IHW09117', 'IHW09118', 'IHW09122', 'IHW09125', 'IHW09175', 'IHW09198', 'IHW09200', 'IHW09224', 'IHW09251', 'IHW09359', 'IHW09409']

# Original 16 samples
hprc_samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]
ihw_samples = ['IHW09021', 'IHW09049', 'IHW09071', 'IHW09117', 'IHW09118', 'IHW09122', 'IHW09125', 'IHW09175', 'IHW09198', 'IHW09200', 'IHW09224', 'IHW09245', 'IHW09251', 'IHW09359', 'IHW09364', 'IHW09409']

twist_prebuilt_panel = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5']

class_I = ['HLA-A', 'HLA-B', 'HLA-C']
class_II = ['HLA-DPA1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5']
class_III =["RPL15P4", "MCCD1", "DDX39B", "ATP6V1G2-DDX39B", "SNORD117", "SNORD84", "DDX39B-AS1", "ATP6V1G2", "NFKBIL1", "ENSG00000289406", "LTA", "TNF", "LTB", "LST1", "NCR3", "UQCRHP1", "AIF1", "ENSG00000289375", "PRRC2A", "ENSG00000291302", "SNORA38", "ENSG00000289282", "MIR6832", "BAG6", "APOM", "C6orf47", "C6orf47-AS1", "GPANK1", "Y_RNA", "CSNK2B", "ENSG00000263020", "LY6G5B", "LY6G5C", "ABHD16A", "ENSG00000204422", "MIR4646", "LY6G6F", "LY6G6F-LY6G6D", "LY6G6E", "LY6G6D", "MPIG6B", "LY6G6C", "DDAH2", "CLIC1", "MSH5", "MSH5-SAPCD1", "RNU6-850P", "SAPCD1", "SAPCD1-AS1", "VWA7", "VARS1", "Y_RNA", "LSM2", "HSPA1L", "HSPA1A", "ENSG00000289637", "ENSG00000289829", "HSPA1B", "ENSG00000285565", "SNHG32", "SNORD48", "SNORD52", "NEU1", "SLC44A4", "EHMT2-AS1", "EHMT2", "C2", "ZBTB12", "ENSG00000244255", "C2-AS1", "CFB", "NELFE", "MIR1236", "SKIC2", "DXO", "STK19", "C4A", "C4A-AS1", "ENSG00000290788", "CYP21A1P", "TNXA", "STK19B", "C4B", "C4B-AS1", "CYP21A2", "TNXB", "RNA5SP206", "ENSG00000284829", "ENSG00000286974", "ATF6B", "FKBPL", "PRRT1", "ENSG00000285085", "ENSG00000284954", "PPT2", "PPT2-EGFL8", "EGFL8", "AGPAT1", "MIR6721", "RNF5", "MIR6833", "AGER", "ENSG00000273333", "PBX2", "GPSM3", "NOTCH4", "TSBP1-AS1", "TSBP1", "HNRNPA1P2", "RNU6-603P", "BTNL2"]


depth_pass_pacbio = []
prop_20x_pacbio = []
prop_30x_pacbio = []

depth_threshold = 30
prop_threshold = 0.9

# Protein Coding Captured. Boundaries identified from Shiina et al. 2009. 
class_I = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G']


def process_coverage_dict():
    global biotype_dict
    with open(biotype_dict_file, "r") as json_file1:
        biotype_dict = json.load(json_file1)

    with open(input_json_file, "r") as json_file2:
        coverage_dict = json.load(json_file2)

    global pacbio_gene_dict 
    pacbio_gene_dict = coverage_dict['revio']['gene']
    #pacbio_gene_dict = coverage_dict['promethion']['gene']

    for gene, data in pacbio_gene_dict.items():
        depth_list = []
        prop_20x_list = []
        prop_30x_list = []
        for sample in sample_list:
            depth_list.append(data[sample]['coverage_depth'])
            prop_20x_list.append(data[sample]['prop_20x'])
            prop_30x_list.append(data[sample]['prop_30x'])
        
        if min(depth_list) > depth_threshold:
            depth_pass_pacbio.append(gene)
        
        else:
            gene_name = data['name']
            failed_samples = dict()
            
            for sample in sample_list:
                if data[sample]['coverage_depth'] < depth_threshold:
                    failed_samples[sample] = data[sample]['coverage_depth']
            
            print("{} ({}) had {} samples with mean coverage depth below 30X".format(gene_name, gene, len(failed_samples)))
            for key,value in failed_samples.items():
                print(key,value)

        if min(prop_20x_list) >= prop_threshold:
            prop_20x_pacbio.append(gene)
        if min(prop_30x_list) >= prop_threshold:
            prop_30x_pacbio.append(gene)

def compare_groups():
    depth_dict = {
    "all": {"Class_I": [], "Class_II": [], "Class_III": []},
    "hprc": {"Class_I": [], "Class_II": [], "Class_III": []},
    "ihw": {"Class_I": [], "Class_II": [], "Class_III": []}
    }

    for gene in depth_pass_pacbio:
        hprc_depth_values = []
        ihw_depth_values = []
        all_depth_values = []
        for sample in hprc_samples:
            hprc_depth_values.append(pacbio_gene_dict[gene][sample]['coverage_depth'])
            all_depth_values.append(pacbio_gene_dict[gene][sample]['coverage_depth'])
        for sample in ihw_samples:
            ihw_depth_values.append(pacbio_gene_dict[gene][sample]['coverage_depth'])
            all_depth_values.append(pacbio_gene_dict[gene][sample]['coverage_depth'])
        
        mean_depth_hprc = mean(hprc_depth_values)
        mean_depth_ihw = mean(ihw_depth_values)
        mean_depth_all = mean(all_depth_values)
       
        if pacbio_gene_dict[gene]['name'] in class_I:
            depth_dict["hprc"]["Class_I"].append(mean_depth_hprc)
            depth_dict["ihw"]["Class_I"].append(mean_depth_ihw)
            depth_dict["all"]["Class_I"].append(mean_depth_all)
        elif pacbio_gene_dict[gene]['name'] in class_II:
            depth_dict["hprc"]["Class_II"].append(mean_depth_hprc)
            depth_dict["ihw"]["Class_II"].append(mean_depth_ihw)
            depth_dict["all"]["Class_II"].append(mean_depth_all)
        elif pacbio_gene_dict[gene]['name'] in class_III:
            depth_dict["hprc"]["Class_III"].append(mean_depth_hprc)
            depth_dict["ihw"]["Class_III"].append(mean_depth_ihw)
            depth_dict["all"]["Class_III"].append(mean_depth_all)

    hprc_c1 = mean(depth_dict["hprc"]["Class_I"])
    hprc_c2 = mean(depth_dict["hprc"]["Class_II"])
    hprc_c3 = mean(depth_dict["hprc"]["Class_III"])

    ihw_c1 = mean(depth_dict["ihw"]["Class_I"])
    ihw_c2 = mean(depth_dict["ihw"]["Class_II"])
    ihw_c3 = mean(depth_dict["ihw"]["Class_III"])

    all_c1 = mean(depth_dict["all"]["Class_I"])
    all_c2 = mean(depth_dict["all"]["Class_II"])
    all_c3 = mean(depth_dict["all"]["Class_III"])

    print("Mean HPRC Class I coverage depth: {:.0f}".format(hprc_c1))
    print("Mean HPRC Class II coverage depth: {:.0f}".format(hprc_c2))
    print("Mean HPRC Class III coverage depth: {:.0f}".format(hprc_c3))

    print("Mean IHW Class I coverage depth: {:.0f}".format(ihw_c1))
    print("Mean IHW Class II coverage depth: {:.0f}".format(ihw_c2))
    print("Mean IHW Class III coverage depth: {:.0f}".format(ihw_c3))

    print("HPRC Enrichment of Class I over class III: {:.1f}X".format(hprc_c1/hprc_c3))
    print("HPRC Enrichment of Class II over Class III: {:.1f}X".format(hprc_c2/hprc_c3))
    print("HPRC Enrichment of Class I/II over Class III: {:.1f}".format(((hprc_c1+hprc_c2)/2)/hprc_c3))

    print("IHW Enrichment of Class I over class III: {:.1f}X".format(ihw_c1/ihw_c3))
    print("IHW Enrichment of Class II over Class III: {:.1f}X".format(ihw_c2/ihw_c3))
    print("IHW Enrichment of Class I/II over Class III: {:.1f}".format(((ihw_c1+ihw_c2)/2)/ihw_c3))

    print("IHW Enrichment over HPRC (Class I): {:.1f}".format(ihw_c1/hprc_c1))
    print("IHW Enrichment over HPRC (Class II): {:.1f}".format(ihw_c2/hprc_c2))
    print("IHW Enrichment over HPRC (Class III): {:.1f}".format(ihw_c3/hprc_c3))

    print("All Samples Class I coverage depth: {:.1f}".format(all_c1))
    print("All Samples Class II coverage depth: {:.1f}".format(all_c2))
    print("All Samples Class III coverage depth: {:.1f}".format(all_c3))

    print("All Samples Enrichment of Class I over class III: {:.1f}X".format(all_c1/all_c3))
    print("All Samples Enrichment of Class II over Class III: {:.1f}X".format(all_c2/all_c3))
    print("All Samples Enrichment of Class I/II over Class III: {:.1f}".format(((all_c1+all_c2)/2)/all_c3))

def print_output():
    print("{} genes had a mean depth of 30X or more across all samples".format(len(depth_pass_pacbio)))
    print("{} genes had all of their bases covered by at least 20 reads".format(len(prop_20x_pacbio)))
    print("{} genes had all of their bases covered by at least 30 reads".format(len(prop_30x_pacbio)))

    targeted_genes = set([gene for item in open(targeted_genes_file, "r").read().splitlines() for gene in item.split("\t")[3].split(";")])
    targeted_genes.update(twist_prebuilt_panel)
    # Ignore HLA-DRB3 and HLA-DRB4 for now because they're not in the hg 38 gtf/gff
    targeted_genes.difference_update({"HLA-DRB3", "HLA-DRB4"})

    # Define captured genes as _X depth across all samples
    captured_genes_depth = []
    for gene in depth_pass_pacbio:
        gene_name = pacbio_gene_dict[gene]["name"]
        captured_genes_depth.append(gene_name)

    captured_genes_20x = []
    for gene in prop_20x_pacbio:
        gene_name = pacbio_gene_dict[gene]["name"]
        captured_genes_20x.append(gene_name)

    if targeted_genes.issubset(captured_genes_depth):
        print("All {} targeted genes had a mean coverage depth greater than {}X across all samples".format(len(targeted_genes), depth_threshold))
    else:
        overlap = len(set(captured_genes_depth).intersection(targeted_genes))
        failed = set(targeted_genes).difference(captured_genes_depth)
        print("{}/{} targeted genes had a mean coverage depth greater than {}X across all samples".format(overlap, len(targeted_genes), depth_threshold))
        print("{} failed genes: {}".format(len(failed),failed))

    if targeted_genes.issubset(captured_genes_20x):
        print("All {} targeted genes had all their bases covered by 20 reads across all samples".format(len(targeted_genes)))
    else:
        overlap = len(set(captured_genes_20x).intersection(targeted_genes))
        failed = set(targeted_genes).difference(captured_genes_20x)
        print("{}/{} targeted genes did not have all bases covered by 20 reads across all samples".format(overlap, len(targeted_genes)))
        print("{} failed genes: {}".format(len(failed),failed))


    captured_protein_coding = 0
    for gene in depth_pass_pacbio:
        if biotype_dict[gene] == "protein_coding":
            captured_protein_coding += 1
        else:
            name = pacbio_gene_dict[gene]["name"]
            #print(name)
    #print(captured_protein_coding)

def main():
    process_coverage_dict()
    compare_groups()
    print_output()

if __name__ == "__main__":
    main()