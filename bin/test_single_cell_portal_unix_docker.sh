today=`date | cut -d " " -f 3,2,4 --output-delimiter=_ | cut -d ":" -f 1,2,3 --output-delimiter=_`
resources=test_data/resources
run_output=testing_pipeline_${today}
samples=testing_data/samples

python scripts/single_cell_analysis.py \
    --gtf ${resources}/hg19_22_abridged.gtf \
    --log ${run_output}/test_pipeline.log \
    --no_trim \
    --out_dir ${run_output} \
    --reference_fasta ${resources}/hg19_22_abridged.fa \
    --reference_flat ${resources}/hg19_22_abridged.refFlat \
    --rrna_intervals ${resources}/hg19_22_abridged.rRNA.intervals \
    --rsem_index ${resources}/hg19_22_abridged \
    --sample_left ${samples}/no_barcodes_left.fq \
    --sample_right ${samples}/no_barcodes_right.fq \
    --update "trimmomatic.jar:/seq/regev_genome_portal/SOFTWARE/Trimmomatic"
