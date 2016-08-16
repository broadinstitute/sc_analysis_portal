workflow scone {

    Int in_threads
    File in_reference
    File in_gtf
    File in_left_sample
    File in_ref_flat
    File in_ribosome_interval
    File in_right_sample
    File in_index
    String in_library
    String in_machine
    String in_output_dir
    String in_picard_jars
    String in_plateform
    String in_sample_name

    call rsem {
        input: dir=${in_output_dir},
               name=${in_sample_name},
               threads=${in_threads},
               left_sample=${in_left_sample},
               right_sample=${in_right_sample},
               index=${in_index}
    }

    call add_rg {
        input: dir=${in_output_dir},
               library=${in_library},
               machine=${in_machine},
               name=${in_sample_name},
               path=${in_picard_jars},
               plateform=${in_plateform},
               bam=${rsem.out_bam}
    }

    call reorder_bam {
        input: dir=${in_output_dir},
               input_bam=${gatk_add_rg.out_bam},
               path=${in_picard_jars},
               reference=${in_reference},
               name=${in_sample_name}
    }

    call mark_duplicates {
        input: dir=${in_output_dir},
               input_bam=${reorder_bam.out_bam},
               name=${in_sample_name},
               path=${in_picard_jars},
               reference=${in_reference}
    }

    call rnaseqqc {
        input: dir=${in_output_dir},
               input_bam=${mark_duplicates.out_bam},
               gtf=${in_gtf},
               name=${in_sample_name},
               path=${in_rnaseqqc_jars},
               reference=${in_reference}
    }

    call fastqc {
        input: dir=${in_output_dir},
               input_bam=${mark_duplicates.out_bam},
               gtf=${in_gtf},
               name=${in_sample_name},
               path=${in_fastqc_jars},
               reference=${in_reference}
    }

    call collect_rna_seq_metrics {
        input: dir=${in_output_dir}/qc,
               input_bam=${mark_duplicates.out_bam},
               path=${in_picard_jars},
               ref_flat=${in_ref_flat},
               ribosome_intervals=${in_ribosome_intervals},
               name=${in_sample_name}
    }

    call estimate_library_complexity {
        input: dir=${in_output_dir}/qc,
               input_bam=${mark_duplicates.out_bam},
               path=${in_picard_jars},
               name=${in_sample_name}
    }

    call collect_alignment_summary_metrics {
        input: dir=${in_output_dir}/qc,
               input_bam=${mark_duplicates.out_bam},
               path=${in_picard_jars},
               reference=${in_reference},
               name=${in_sample_name}
    }

    call collect_insert_size_metrics {
        input: dir=${in_output_dir}/qc,
               input_bam=${mark_duplicates.out_bam},
               path=${in_picard_jars},
               name=${in_sample_name}
    }
}

task add_rg {
    String library
    String machine
    String name
    String path
    String plateform
    File bam
    command {
        java -jar ${path}/AddOrReplaceReadGroups.jar \
            I=${bam} \
            O=${out_bam} \
            SO=coordinate \
            RGLB=${library} \
            RGPL=${plateform} \
            RGPU=${machine} \
            RGSM=${name}
    }
    output {
        File out_bam = "${name}_rg.bam"
    }
}

task mark_duplicates {
    String path
    String name
    File in_bam
    command {
        java -jar ${path}/MarkDuplicates.jar \
            I=${in_bam} \
            O=${out_bam} \
            CREATE_INDEX=true \
            M=${qc_metrics}
    }
    output {
        File out_bam = "${name}_dedubbed.bam"
        File qc_metrics = "${name}_qc_metrics.txt"
    }
}

task reorder_bam {
    File bam
    File reference
    String path
    String name
    command {
        java -Xmx4G -jar ${path}/ReorderSam.jar \
            I= ${bam}\
            R= ${reference}\
            O= ${out_bam}
    }
    output {
        File out_bam = "${name}_reorder.bam"
    }
}

task rsem {
    Int threads
    File dir
    File left_sample
    File right_sample
    File index
    String name
    command {
        rsem-calculate-expression \
            -p ${threads} \
            --paired-end \
            --bowtie2 \
            --estimate-rspd \
            --output-genome-bam \
            ${left_sample} ${right_sample} \
            ${index} ${dir}/${name}/RSEM
    }
    output {
        File out_gene_results = ${dir}/${name}/RSEM.genes.results
        File out_isoform_results = ${dir}/${name}/RSEM.isoform.results
        File out_bam = ${dir}/${name}/RSEM.genome.bam
        File out_sorted_bam = ${dir}/${name}/RSEM.genome.sorted.bam
        File out_sorted_index = ${dir}/${name}/RSEM.genome.sorted.bai
        File out_transcript_bam = ${dir}/${name}/RSEM.transcript.bam
        File out_trans_sorted_bam = ${dir}/${name}/RSEM.transcript.sorted.bam
        File out_trans_sorted_bai = ${dir}/${name}/RSEM.transcript.sorted.bam.bai
        File out_model_count = ${dir}/${name}/RSEM.stat/RSEM.cnt
        File out_model = ${dir}/${name}/RSEM.stat/RSEM.model
        File out_model_theta = ${dir}/${name}/RSEM.stat/RSEM.theta
    }
}

task rnaseqqc {
    File input_bam
    File dir
    File reference
    File gtf
    String path
    String name
    command {
        java -Xmx4g -jar ${path}/RNA-SeQC.jar \
            -n 1000 \
            -s "sample|"${bam}"|RNASeQC" \
            -t ${gtf} \
            -r ${fasta} \
            -o $${dir}/{name}/rnaseqqc
    }
    output {
        File out_genome_qc = "${name}/RSEM_genome_qc.txt"
        File out_isoform_results = "${name}/RSEM.isoforms.results"
        File out_gene_results = "${name}/RSEM.genes.results"
        File out_sample_metrics = "${name}/rnaseqqc/metrics.tsv"
    }

task fastqc {

    String dir
    String name
    File fastq_left
    File fastq_right

    command {
        fastqc ${fastq_left} ${fastq_right} \
            --outdir=${dir}/${name}
    }
    output {
        File out=${dir}/${name}
    }
}

task collect_rna_seq_metrics {

    String dir
    String name
    String path
    File input_bam
    File ref_flat
    File ribosome_intervals

    command {
        java -jar ${path}/CollectRnaSeqMetrics.jar \
            I=${input_bam} \
            O=${dir}/${name}/collect_rna_seq_metrics.txt \
            REF_FLAT=${ref_flat} \
            STRAND=NONE \
            RIBOSOMAL_INTERVALS=${ribosome_intervals}
    }
    output {
        File out = ${dir}/${name}/collect_rna_seq_metrics.txt
    }
}

task estimate_library_complexity {

    String dir
    String name
    String path
    File input_bam

    command {
        java -jar ${dir}/EstimateLibraryComplexity.jar \
            I=${input_bam} \
            O=${dir}/${name}/estimate_library_complexity.txt
    }
    output {
        File out = ${dir}/${name}/estimate_library_complexity.txt
    }
}

task collect_insert_size_metrics {

    String dir
    String name
    String path
    File input_bam
    command {
        java -jar ${path}/CollectInsertSizeMetrics.jar \
            I=${input_bam} \
            O=${dir}/${name}/collect_insert_size.txt \
            H=${dir}/${name}/collect_insert_size.pdf \
            M=0.5
    }
    output {
        File out = ${dir}/${name}/collect_insert_size.txt
        File out_pdf = ${dir}/${name}/collect_insert_size.pdf 
    }
}

task collect_alignment_summary_metrics {

    String dir
    String name
    String path
    File input_bam
    File reference
    command {
        java -jar ${path}/CollectAlignmentSummaryMetrics.jar \
        REFERENCE_SEQUENCE=${reference} \
        INPUT=${input_bam} \
        OUTPUT=${dir}/${name}/collect_aligment_summary_metrics.txt
    }
    output {
        File out = ${dir}/${name}/collect_aligment_summary_metrics.txt
    }
}
