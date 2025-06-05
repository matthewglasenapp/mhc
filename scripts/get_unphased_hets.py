import pysam

vcf_path = "/hb/groups/cornejo_lab/matt/hla_capture/pacbio/phased_vcf_hiphase/HG002.dedup.trimmed.hg38.chr6.phased.joint.vcf.gz"
output_bed = "unphased_hets.bed"

vcf = pysam.VariantFile(vcf_path)
region = "chr6"
start = 28000000
end = 34000000

with open(output_bed, "w") as out:
    for record in vcf.fetch(region, start, end):
        for sample in record.samples:
            gt = record.samples[sample]['GT']
            phased = record.samples[sample].phased
            if gt and len(gt) == 2 and gt[0] != gt[1] and not phased:
                bed_start = record.pos - 1
                bed_end = record.pos
                out.write(f"{record.chrom}\t{bed_start}\t{bed_end}\n")
