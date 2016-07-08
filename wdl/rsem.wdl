workflow scone {
    Int in_threads
    File in_reference
    File in_gtf
    File in_left_sample
    File in_right_sample
    File in_index
    String in_library
    String in_machine
    String in_plateform
    String in_sample_name

    call rsem {
        input: name=${in_sample_name}, threads=${in_threads}, left_sample=${in_left_sample}, right_sample=${in_right_sample}, index=${in_index}
    }

    call gatk_add_rg {
        input: library=${in_library}, machine=${in_machine}, name=${in_sample_name}, plateform=${in_plateform}, bam=${rsem.out_bam}
    }

    call reorder_bam {
        input: bam=${gatk_add_rg.out_bam}, reference=${in_reference}, name=${in_sample_name}
    }

    call mark_duplicates {
        input bam=${reorder_bam.out_bam}, reference=${in_reference}, name=${in_sample_name}
    }

    call rnaseqqc {
        input: bam=${mark_duplicates.out_bam}, name=${in_sample_name}, reference=${in_reference}, gtf=${in_gtf}
 
   }
}

task gatk_add_rg {
    String library
    String machine
    String name
    String plateform
    File bam
    command {
        java -jar AddOrReplaceReadGroups.jar \
            I=${bam} \
            O=${out_bam} \
            SO=coordinate \
            RGLB=${library} \
            RGPL=${plateform} \
            RGPU=${machine} \
            RGSM=${name}
    }
    output {
        File out_bam = ${name}_rg.bam
    }
}

task mark_duplicates {
    String name
    File in_bam
    command {
        java -jar MarkDuplicates.jar \
            I=${in_bam} \
            O=${out_bam} \
            CREATE_INDEX=true \
            M=${qc_metrics}
    }
    output {
        File out_bam = ${name}_dedubbed.bam 
        File qc_metrics = ${name}_qc_metrics.txt
    }
}

task reorder_bam {
    File bam
    File reference
    String name
    command {
        java -Xmx4G -jar ReorderSam.jar \
            I= ${bam}\
            R= ${reference}\
            O= ${out_bam}
    }
    output {
        File out_bam = ${name}_reorder.bam
    }
}

task rsem {
    Int threads
    File left_sample
    File right_sample
    File index
    String name
    command {
        rsem-calculate-expression \
            --bowtie-chunkmbs 512 \
            -p ${threads} \
            --paired ${left_sample} ${right_sample} \
            ${index} ${name}/RSEM
    }
    output {
        File out_gene_results = ${name}/RSEM.genes.results
        File out_isoform_results = ${name}/RSEM.isoform.results
        File out_bam = ${name}/RSEM.genome.bam
        File out_sorted_bam = ${name}/RSEM.genome.sorted.bam
        File out_sorted_index = ${name}/RSEM.genome.sorted.bai
        File out_transcript_bam = ${name}/RSEM.transcript.bam
        File out_trans_sorted_bam = ${name}/RSEM.transcript.sorted.bam
        File out_trans_sorted_bai = ${name}/RSEM.transcript.sorted.bam.bai
        File out_model_count = ${name}/RSEM.stat/RSEM.cnt
        File out_model = ${name}/RSEM.stat/RSEM.model
        File out_model_theta = ${name}/RSEM.stat/RSEM.theta
    }
}
task rsem {
    Int threads
    File left_sample
    File right_sample
    File index
    String name
    command {
        rsem-calculate-expression \
            -p ${threads} \
            --paired-end \
            --bowtie2 \
            --bowtie2-path ${bowtie2path} \
            --estimate-rspd \
            --fragment-length-max \ # Paired 
            --output-genome-bam \
            ${left_sample} ${right_sample} \
            ${index} ${name}/RSEM
    }
    output {
        File out_gene_results = ${name}/RSEM.genes.results
        File out_isoform_results = ${name}/RSEM.isoform.results
        File out_bam = ${name}/RSEM.genome.bam
        File out_sorted_bam = ${name}/RSEM.genome.sorted.bam
        File out_sorted_index = ${name}/RSEM.genome.sorted.bai
        File out_transcript_bam = ${name}/RSEM.transcript.bam
        File out_trans_sorted_bam = ${name}/RSEM.transcript.sorted.bam
        File out_trans_sorted_bai = ${name}/RSEM.transcript.sorted.bam.bai
        File out_model_count = ${name}/RSEM.stat/RSEM.cnt
        File out_model = ${name}/RSEM.stat/RSEM.model
        File out_model_theta = ${name}/RSEM.stat/RSEM.theta
    }
}

task rnaseqqc {
    File bam
    File reference
    File gtf
    String name
    command {
        java -Xmx4g -jar RNA-SeQC.jar \
            -n 1000 \
            -s "sample|"${bam}"|RNASeQC" \
            -t ${gtf} \
            -r ${fasta} \
            -o ${name}/rnaseqqc
    }
    output {
        File out_qc = 
    }

