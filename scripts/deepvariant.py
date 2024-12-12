import os

wd = "/hb/scratch/ogarci12/deepvariant"
# Revio data with duplicates removed
input_dir = "/hb/scratch/ogarci12/PacBio_mm2/rm_dup/"

# Revio Output Dir
output_dir = "/hb/scratch/ogarci12/deepvariant/Revio_VCF/"

# PromethION BAM with Duplicates Removed
promethion_input_dir = "/hb/scratch/mglasena/deepvariant_whathap/Promethion_mm2/mapped_raw/"

# PromethION/Clair3 output dir
clair3_output_dir = "/hb/scratch/mglasena/deepvariant_whathap/PromethION_VCF/"

# Output dir for whatshap
phased_bam_dir = "/hb/scratch/ogarci12/deepvariant/Revio_Phased_BAM"

# Reference fasta
ref = "/hb/scratch/ogarci12/deepvariant/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

deepvariant_sif = "/hb/scratch/ogarci12/deepvariant/deepvariant.sif"
clair3_sif = "/hb/scratch/mglasena/deepvariant_whathap/clair3_latest.sif"
clair3_model_path = "/hb/home/mglasena/.conda/envs/clair3/bin/models/r941_prom_sup_g5014"

# Threads for Clair3
threads = 8

make_output_dir = "mkdir -p {}".format(output_dir)
os.system(make_output_dir)

make_phased_bam_dir = "mkdir -p {}".format(phased_bam_dir)
os.system(make_phased_bam_dir)

samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01891", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", "IHW09409", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]

def index_bam(sample):
    input_file = input_dir + sample + "_rm_duplicates.bam"
    index = "samtools index {}".format(input_file)
    os.system(index)

def run_deepvariant(sample):
    os.chdir(output_dir)
    run_deepvariant = "singularity exec --bind {}/:/data --bind {}/:/input --bind /hb/scratch/ogarci12/deepvariant:/reference {} /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref=/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa --reads=/input/{}_rm_duplicates.bam --output_vcf=/data/{}.hg38.Pacbio.vcf.gz --output_gvcf=/data/{}.hg38.Pacbio.g.vcf.gz --num_shards=8".format(output_dir, input_dir, deepvariant_sif, sample, sample, sample)
    os.system(run_deepvariant)

def run_clair3(sample):
    input_file = promethion_input_dir + sample + ".hg38.minimap.promethion.raw.bam"
    print("Variant Calling {}".format(input_file))
    index = "samtools index {}".format(input_file)
    os.system(index)

    run_clair3 = "run_clair3.sh --bam_fn={} --ref_fn={} --platform=ont --model_path={} --output={} --threads={}".format(input_file, ref, clair3_model_path, clair3_output_dir, threads)
    os.system(run_clair3)

def index_vcf(sample):
    os.chdir(output_dir)
    tabix = "tabix {}.hg38.Pacbio.vcf.gz".format(sample)
    index = "bcftools index {}.hg38.Pacbio.vcf.gz".format(sample)
    os.system(tabix)
    os.system(index)

def run_whatshap(sample):
    run_whatshap = "whatshap phase --output {}/{}.hg38.Pacbio.phased.vcf.gz --reference={} {}/{}.hg38.Pacbio.vcf.gz {}/{}.pacbio.primary.hg38.bam".format(output_dir, sample, ref, output_dir, sample, input_dir, sample)
    os.system(run_whatshap)
    index = "bcftools index {}/{}.hg38.Pacbio.phased.vcf.gz".format(output_dir, sample)
    tabix = "tabix {}/{}.hg38.Pacbio.phased.vcf.gz".format(output_dir, sample)
    os.system(index)
    os.ystem(tabix)

    run_haplotag = "whatshap haplotag -o {}/{}.hg38.Pacbio.phased.bam --reference={} {}/{}.hg38.Pacbio.vcf.gz {}/{}.pacbio.primary.hg38.bam".format(phased_bam_dir, sample, ref, output_dir, sample, input_dir, sample)
    os.system(run_haplotag)

def change_permissions():
    os.system("chmod -R 777 {}".format(wd))
    os.system("chmod -R 777 {}".format(input_dir))

def main():
    array_id = os.environ["array_id"]
    print("Array ID: {}".format(array_id))
    sample = samples[int(array_id)]
    print(sample)
    #index_bam(sample)
    #run_deepvariant(sample)
    run_clair3(sample)
    #index_vcf(sample)
    #run_whatshap(sample)
    #change_permissions()

if __name__ == "__main__":
    main()
