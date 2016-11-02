today=`date | cut -d " " -f 3,2,4 --output-delimiter=_ | cut -d ":" -f 1,2,3 --output-delimiter=_`
mkdir testing_prep_${today}

python scripts/prep_resources.py \
    --reference test_data/resources/hg19_22_abridged.fa \
    --gtf test_data/resources/hg19_22_abridged.gtf \
    --out_dir testing_prep_${today} \
    --log testing_prep_${today}/test_prep.log \
    --update "picard.jar:/seq/regev_genome_portal/SOFTWARE/Picard/picard-tools-2.0.1"
