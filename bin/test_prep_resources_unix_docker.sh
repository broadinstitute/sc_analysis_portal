IMAGE=4c78f1d77a87
resources=/home/ubuntu/data
output_dir=/home/ubuntu/sc_analysis_portal/test_prep_run

sudo docker run \
    -v ${resources}:/usr/local/data \
    -v ${output_dir}:/usr/local/data/test \
    ${IMAGE} python sc_analysis_portal/scripts/prep_resources.py \
    --reference_fasta /usr/local/data/hg19.fasta \
    --gtf /usr/local/data/hg19.gtf \
    --out_dir /usr/local/data/test \
    --log /usr/local/data/test/test_run.log \
    --update "picard.jar:/usr/local/bin"
