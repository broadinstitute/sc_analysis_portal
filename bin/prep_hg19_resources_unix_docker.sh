python scripts/prep_resources.py \
    --reference /ahg/regevdata/projects/scanalysisportal/data/from_dropseq/hg19.fa \
    --gtf /ahg/regevdata/projects/scanalysisportal/data/from_dropseq/hg19.gtf \
    --out_dir prep_hg19 \
    --log prep_hg19/prep.log \
    --update "picard.jar:/seq/regev_genome_portal/SOFTWARE/Picard/picard-tools-2.0.1"
