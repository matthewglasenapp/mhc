import os

# Revio mapped BAM with duplicates removed
revio_input_dir = "/hb/scratch/mglasena/MHC/mapped_bam/revio/"

# Revio Output Dir
revio_output_dir = "/hb/scratch/mglasena/MHC/genotypes/revio/"

# PromethION mapped BAM with duplicates removed
promethion_input_dir = "/hb/scratch/mglasena/MHC/mapped_bam/promethion/"

# PromethION output dir
promethion_output_dir = "/hb/scratch/mglasena/MHC/genotypes/promethion/"

# Reference fasta
ref = "/hb/scratch/mglasena/MHC/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

deepvariant_sif = "/hb/scratch/mglasena/MHC/deepvariant_sif/deepvariant.sif"
clair3_model_path = "/hb/home/mglasena/.conda/envs/clair3/bin/models/r941_prom_sup_g5014"

# Threads for Clair3
threads = 8

make_revio_output_dir = "mkdir -p {}".format(revio_output_dir)
os.system(make_revio_output_dir)

make_promethion_output_dir = "mkdir -p {}".format(promethion_output_dir)
os.system(make_promethion_output_dir)

samples = ["HG002", "HG003", "HG004", "HG005", "HG01106", "HG01258", "HG01891", "HG01928", "HG02055", "HG02630", "HG03492", "HG03579", "IHW09021", "IHW09049", "IHW09071", "IHW09117", "IHW09118", "IHW09122", "IHW09125", "IHW09175", "IHW09198", "IHW09200", "IHW09224", "IHW09245", "IHW09251", "IHW09359", "IHW09364", "IHW09409", "NA19240", "NA20129", "NA21309", "NA24694", "NA24695"]

def run_deepvariant(sample):
    os.chdir(revio_output_dir)
    run_deepvariant = "singularity exec --bind {}/:/data --bind {}/:/input --bind /hb/scratch/mglasena/MHC/reference:/reference {} /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref=/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa --reads=/input/{}.hg38.revio_rmdup.bam --output_vcf=/data/{}.hg38.revio.vcf.gz --output_gvcf=/data/{}.hg38.revio.g.vcf.gz --num_shards=8".format(revio_output_dir, revio_input_dir, deepvariant_sif, sample, sample, sample)
    os.system(run_deepvariant)

def run_clair3(sample):
    input_file = promethion_input_dir + sample + ".hg38.promethion_rmdup.bam"
    output_dir = promethion_output_dir + sample
    os.system("mkdir -p " + output_dir)
    print("Variant Calling {}".format(input_file))

    index = "samtools index {}".format(input_file)
    os.system(index)

    run_clair3 = "run_clair3.sh --bam_fn={} --ref_fn={} --platform=ont --model_path={} --output={} --threads={}".format(input_file, ref, clair3_model_path, output_dir, threads)
    os.system(run_clair3)

def index_vcf(sample):
    os.chdir(revio_output_dir)
    tabix = "tabix {}.hg38.revio.vcf.gz".format(sample)
    index = "bcftools index {}.hg38.revio.vcf.gz".format(sample)
    os.system(tabix)
    os.system(index)

def main():
    array_id = os.environ["array_id"]
    print("Array ID: {}".format(array_id))
    sample = samples[int(array_id)]
    print(sample)
    #run_deepvariant(sample)
    run_clair3(sample)
    #index_vcf(sample)

if __name__ == "__main__":
    main()
