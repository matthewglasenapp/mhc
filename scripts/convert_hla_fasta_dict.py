import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

#input_json = "/Users/matt/Desktop/fasta_dict.json"
#input_json = "/Users/matt/Desktop/cds_dict.json"
output_fasta = "/Users/matt/Desktop/HLA_Class_I_haplotypes.fa"

def json_to_fasta(input_json, output_fasta):
    records = []

    with open(input_json) as f:
        fasta_data = json.load(f)

    for platform, genes in fasta_data.items():
        for gene, samples in genes.items():
            for sample, haplotypes in samples.items():
                hap1_name = f"{sample}_{gene}_{platform}_1"
                hap2_name = f"{sample}_{gene}_{platform}_2"
                hap1_seq = haplotypes[0]
                hap2_seq = haplotypes[1]

                records.append(SeqRecord(Seq(hap1_seq), id=hap1_name, description=""))
                records.append(SeqRecord(Seq(hap2_seq), id=hap2_name, description=""))

    SeqIO.write(records, output_fasta, "fasta")
    print(f"Wrote {len(records)} records to {output_fasta}")

if __name__ == "__main__":
    json_to_fasta(input_json, output_fasta)
