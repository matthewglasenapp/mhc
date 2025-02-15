module load dorado/0.9.1

basecalled_reads="/hb/groups/cornejo_lab/HLA_hybrid_capture/06_25_24_R1041_LIG_Cornejo_EXP26/Cornejo/06_25_24_R1041_LIG_Cornejo_EXP26_1_drd0.7.2_sup5.0.0.bam"
outdir="/hb/scratch/mglasena/test_ont/"

dorado demux -o $outdir \
    --emit-summary \
    --barcode-arrangement /hb/scratch/mglasena/test_ont/sample_barcode_arrangement.txt \
    --barcode-sequences /hb/scratch/mglasena/test_ont/sample_barcodes.fa \
    --kit-name MY_CUSTOM_KIT \
    --threads 16 \
    --barcode-both-ends \
    $basecalled_reads



