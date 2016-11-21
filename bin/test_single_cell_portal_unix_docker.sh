IMAGE=4c78f1d77a87
resources=/home/ubuntu/data_for_run
run_output=/home/ubuntu/testing_pipeline_run
samples=/home/ubuntu/data

sudo docker run \
    -v ${resources}:/usr/local/data \
    -v ${run_output}:/usr/local/data/test \
    -v ${samples}:/usr/local/data/samples \
    ${IMAGE} single_cell_analysis.py \
    --gtf /usr/local/data/hg19.gtf \
    --no_trim \
    --out_dir /usr/local/data/test \
    --reference_fasta /usr/local/data/hg19.fasta \
    --reference_flat /usr/local/data/hg19.refFlat \
    --rrna_intervals /usr/local/data/hg19.rRNA.intervals \
    --rsem_index /usr/local/data/hg19 \
    --sample_left /usr/local/data/samples/Left_10000.fq \
    --sample_right /usr/local/data/samples/Right_10000.fq \
    --update "trimmomatic.jar:/usr/local/bin/Trimmomatic" \
    --test

sudo docker run -it \
    -v ${resources}:/usr/local/data \
    -v ${run_output}:/usr/local/data/test \
    -v ${samples}:/usr/local/data/samples \
    ${IMAGE} bash

# --log /usr/local/data/test/test_pipeline_run.log \
