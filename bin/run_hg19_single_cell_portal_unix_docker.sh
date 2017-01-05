resources=/ahg/regevdata/projects/scanalysisportal/data/prep_hg19_22
run_output=run_hg19
samples=/seq/regev_genome_portal/KCO_TESTING/general/rnaseq/breast_cancer

python scripts/single_cell_analysis.py \
    --gtf ${resources}/hg19_22.gtf \
    --log ${run_output}/run_hg19_22.log \
    --no_trim \
    --out_dir ${run_output} \
    --reference_fasta ${resources}/hg19_22.fa \
    --reference_flat ${resources}/hg19_22.refFlat \
    --rrna_intervals ${resources}/hg19_22.rRNA.intervals \
    --rsem_index ${resources}/hg19_22 \
    --sample_right ${samples}/MCF-7/Right/fq.gz \
    --sample_left ${samples}/MCF-7/Left.fq.gz \
    --update "trimmomatic.jar:/seq/regev_genome_portal/SOFTWARE/Trimmomatic,AddOrReplaceReadGroups.jar:/seq/regev_genome_portal/SOFTWARE/Picard/current,ReorderSam.jar:/seq/regev_genome_portal/SOFTWARE/Picard/current,MarkDuplicates.jar:/seq/regev_genome_portal/SOFTWARE/Picard/current,RNA-SeQC.jar:/seq/regev_genome_portal/SOFTWARE/RNASEQ-QC,CollectRnaSeqMetrics.jar:/seq/regev_genome_portal/SOFTWARE/Picard/current,EstimateLibraryComplexity.jar:/seq/regev_genome_portal/SOFTWARE/Picard/current,CollectAlignmentSummaryMetrics.jar:/seq/regev_genome_portal/SOFTWARE/Picard/current,CollectInsertSizeMetrics.jar:/seq/regev_genome_portal/SOFTWARE/Picard/current"
