output_tsv="test.tsv"
output_bed="test.bed"
#input_truth="/hb/groups/cornejo_lab/matt/hla_capture/GIAB_benchmark/HG002/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.vcf.gz"
#input_truth="cleaned_truth.chr6.vcf.gz"
input_truth="/hb/groups/cornejo_lab/matt/hla_capture/GIAB_benchmark/HG002/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.vcf.gz"
#input_query="/hb/groups/cornejo_lab/matt/hla_capture/pacbio/phased_vcf_hiphase/HG002.dedup.trimmed.hg38.chr6.phased.joint.vcf.gz"
input_query="query_cleaned.vcf.gz"

whatshap compare \
--ignore-sample-name \
--tsv-pairwise $output_tsv \
--switch-error-bed $output_bed \
$input_truth \
$input_query

bcftools view -g het -m2 -M2 -v snps -i 'GT="het" && GT~"^[01]\|[01]$"' $input_query -Oz -o query_cleaned.vcf.gz
tabix -p vcf query_cleaned.vcf.gz
