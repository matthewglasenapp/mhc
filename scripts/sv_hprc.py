import subprocess
import shutil
import sys
import os

#==============================
# Download necessary files for analysis

# Download hg19 reference with decoys
# mkdir ref
# Download HG002 paternal haplotype
# wget -P ref https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_021950905.1/ensembl/genome/Homo_sapiens-GCA_021950905.1-hardmasked.fa.gz
# Download HG002 maternal haplotype
# wget -P ref https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_021951015.1/ensembl/genome/Homo_sapiens-GCA_021951015.1-hardmasked.fa.gz
# gunzip ref/Homo_sapiens-GCA_021950905.1-hardmasked.fa.gz
# gunzip ref/Homo_sapiens-GCA_021951015.1-hardmasked.fa.gz

# Map non-ACGT characters to N
# sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' ref/Homo_sapiens-GCA_021950905.1-hardmasked.fa

# sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' ref/Homo_sapiens-GCA_021951015.1-hardmasked.fa

# Download tandem repeat bed annotations
# What to do about this??
#==============================

#==============================
# # MCCD1, BTNL2
# #!genome-build-accession GCA_021951015.1
# 6	ensembl	gene	31464310	31465580	.	+	.	ID=gene:ENSG05525003001;Name=MCCD1;biotype=protein_coding;description=mitochondrial coiled-coil domain 1 [Source:HGNC Symbol%3BAcc:HGNC:20668]%3Bparent_gene_display_xref%3DMCCD1;gene_id=ENSG05525003001;version=1
# 6	ensembl	gene	32329210	32342358	.	-	.	ID=gene:ENSG05525002714;Name=BTNL2;biotype=protein_coding;description=butyrophilin like 2 [Source:HGNC Symbol%3BAcc:HGNC:1142]%3Bparent_gene_display_xref%3DBTNL2;gene_id=ENSG05525002714;version=1

# #!genome-build-accession GCA_021950905.1
# 6	ensembl	gene	31477275	31478545	.	+	.	ID=gene:ENSG05520029244;Name=MCCD1;biotype=protein_coding;description=mitochondrial coiled-coil domain 1 [Source:HGNC Symbol%3BAcc:HGNC:20668]%3Bparent_gene_display_xref%3DMCCD1;gene_id=ENSG05520029244;version=1
# 6	ensembl	gene	32342163	32355332	.	-	.	ID=gene:ENSG05520045463;Name=BTNL2;biotype=protein_coding;description=butyrophilin like 2 [Source:HGNC Symbol%3BAcc:HGNC:1142]%3Bparent_gene_display_xref%3DBTNL2;gene_id=ENSG05520045463;version=1

# # HLA-F, HLA-DPB2
# #!genome-build-accession GCA_021951015.1
# 6	ensembl	gene	29666529	29682266	.	+	.	ID=gene:ENSG05525000195;Name=HLA-F;biotype=protein_coding;description=major histocompatibility complex%2C class I%2C F [Source:HGNC Symbol%3BAcc:HGNC:4963]%3Bparent_gene_display_xref%3DHLA-F;gene_id=ENSG05525000195;version=1
# 6	ensembl	pseudogene	33035199	33051831	.	+	.	ID=gene:ENSG05525002569;Name=HLA-DPB2;biotype=transcribed_unprocessed_pseudogene;description=major histocompatibility complex%2C class II%2C DP beta 2 (pseudogene) [Source:HGNC Symbol%3BAcc:HGNC:4941]%3Bparent_gene_display_xref%3DHLA-DPB2;gene_id=ENSG05525002569;version=1

# #!genome-build-accession GCA_021950905.1
# 6	ensembl	gene	29681374	29697143	.	+	.	ID=gene:ENSG05520019370;Name=HLA-F;biotype=protein_coding;description=major histocompatibility complex%2C class I%2C F [Source:HGNC Symbol%3BAcc:HGNC:4963]%3Bparent_gene_display_xref%3DHLA-F;gene_id=ENSG05520019370;version=1
# 6	ensembl	pseudogene	33144559	33160899	.	+	.	ID=gene:ENSG05520003952;Name=HLA-DPB2;biotype=transcribed_unprocessed_pseudogene;description=major histocompatibility complex%2C class II%2C DP beta 2 (pseudogene) [Source:HGNC Symbol%3BAcc:HGNC:4941]%3Bparent_gene_display_xref%3DHLA-DPB2;gene_id=ENSG05520003952;version=1
#==============================

max_threads = 12

reference_genomes = {
	"hg38": "ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa",
	"hg19": "ref/human_hs37d5.fasta",
	"paternal": "ref/Homo_sapiens-GCA_021950905.1-hardmasked.fa",
	"maternal": "ref/Homo_sapiens-GCA_021951015.1-hardmasked.fa"
}

# tandem_repeat_bed = "ref/human_hs37d5.trf.bed"

fastq_file = "/hb/scratch/mglasena/test_pacbio/processed_data/fastq_rmdup_cutadapt/HG002.dedup.trimmed.fastq.gz"
HG002_RG_string = r'"@RG\tID:m84039_240622_113450_s1\tSM:HG002"'

# Ensure all required tools are installed and executable
def check_required_commands():    
	print("Checking the installation status of the required bioinformatics tools!")

	required_commands = [
		"bgzip",
		"pbmm2",
		"pbsv",
		"samtools",
		"tabix",
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

# Map fastq files to hg19 reference genome
def align_to_reference(reference_fasta, reference_accession):
	if reference_accession == "hg38":
		region = "chr6"
	else:
		region = 6

	output_bam = "mapped_bam/HG002.dedup.trimmed.{reference_accession}.bam".format(reference_accession = reference_accession)

	pbmm2_cmd = "pbmm2 align -j {max_threads} {reference_fasta} --sort --log-level INFO --unmapped --bam-index BAI  {input_file} {output_file} --rg {read_group_string}".format(max_threads = max_threads, reference_fasta = reference_fasta, input_file = fastq_file, output_file = output_bam, read_group_string = HG002_RG_string)

	subprocess.run(pbmm2_cmd, shell=True, check=True)

	# Filter reads for chromosome 6
	filtered_bam = "mapped_bam/HG002.dedup.trimmed.{reference_accession}.chr6.bam".format(reference_accession = reference_accession)
	samtools_cmd = "samtools view -@ {max_threads} -b {input_bam} {region} > {output_bam}".format(max_threads = max_threads, input_bam = output_bam, region = region, output_bam = filtered_bam)

	subprocess.run(samtools_cmd, shell=True, check=True)

	index_cmd = "samtools index {input_bam}".format(input_bam = filtered_bam)

	subprocess.run(index_cmd, shell=True, check=True)

def run_pbsv(reference_fasta, reference_accession):
	if reference_accession == "hg38":
		region = "chr6"
	else:
		region = "6"

	# Run pbsv	
	input_bam = "mapped_bam/HG002.dedup.trimmed.{reference_accession}.chr6.bam".format(reference_accession = reference_accession)
	output_svsig = "pbsv_vcf/HG002.dedup.trimmed.{reference_accession}.chr6.svsig.gz".format(reference_accession = reference_accession)
	output_vcf = "pbsv_vcf/HG002.dedup.trimmed.{reference_accession}.chr6.vcf".format(reference_accession = reference_accession)

	# --tandem-repeats
	pbsv_discover_cmd = "pbsv discover --region {region} {input_bam} {output_svsig}".format(region = region, input_bam = input_bam, output_svsig = output_svsig)
	
	subprocess.run(pbsv_discover_cmd, shell=True, check=True)

	index_cmd = "tabix -c '#' -s 3 -b 4 -e 4 {output_svsig}".format(output_svsig = output_svsig)

	subprocess.run(index_cmd, shell=True, check=True)
	
	pbsv_call_cmd = "pbsv call -j {max_threads} --region {region} --hifi {reference_fasta} {input_svsig} {output_vcf}".format(max_threads = max_threads, region = region, reference_fasta = reference_fasta, input_svsig = output_svsig, output_vcf = output_vcf)

	subprocess.run(pbsv_call_cmd, shell=True, check=True)

	bgzip_cmd = "bgzip {input_file}".format(input_file = output_vcf)
	tabix_cmd = "tabix {input_file}".format(input_file = output_vcf + ".gz")
	subprocess.run(bgzip_cmd, shell=True, check=True)
	subprocess.run(tabix_cmd, shell=True, check=True)

def main():
	check_required_commands()
	
	output_dirs = ["mapped_bam/", "pbsv_vcf/"]
	for dir in output_dirs:
		os.makedirs(dir, exist_ok=True)

	for reference,path in reference_genomes.items():
		reference_fasta = path
		if reference == "hg19" or reference == "hg38":
			reference_accession = reference
		else:
			reference_accession = path.split("/")[1].split("-")[1]

		print(f"Aligning to reference: {reference_fasta} ({reference_accession})")
		align_to_reference(reference_fasta, reference_accession)
		run_pbsv(reference_fasta, reference_accession)

if __name__ == "__main__":
	main()

