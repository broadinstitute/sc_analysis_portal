workflow single_cell_analysis_portal {

    File in_gtf
    File in_output_dir
    File in_reference
    File in_ref_flat
    File in_index
    File in_ribosome_interval
    File in_left_sample
    File in_right_sample
    String in_library
    String in_machine
    String in_plateform
    String in_sample_name
    Int in_threads

    call run_pipeline {
        input: gtf=${in_gtf},
               out_dir=${in_output_dir},
               reference_fasta=${in_reference},
               reference_flat=${in_ref_flat},
               rsem_index=${in_index},
               rrna_intervals=${in_ribosome_interval},
               sample_left=${in_left_sample},
               sample_right=${in_right_sample},
               library=${in_library},
               machine=${in_machine},
               name=${in_name},
               plateform=${in_plateform},
               threads=${in_threads},
    }

task run_pipeline {
    File gtf
    File out_dir
    File reference_fasta
    File reference_flat
    File rsem_index
    File rrna_intervals
    File sample_left
    File sample_right
    String library
    String machine
    String name
    String plateform
    Int threads

    command <<<
        single_cell_analysis.py \
            --trim ${trim}
            --reference_fasta ${reference_fasta} \
            --reference_flat ${reference_flat} \
            --rsem_index ${rsem_index} \
            --rib_intervals ${rrna_intervals} \
            --gtf ${gtf} \
            --sample_left ${sample_left} \
            --sample_right ${sample_right} \
            --add_machine ${machine} \
            --add_name ${name} \
            --add_library ${library} \
            --add_plateform ${plateform} \
            --threads ${threads} \
            --out ${out_dir}
    >>>

    output {
        File trimmed_left fastq = "${out_dir}/"
        File trimmed_right fastq = "${out_dir}/"
        File rsem_bam = "${out_dir}/bam/"
        File prepped_bam = "${out_dir}/bam/"
    }

    runtime {
        docker: "single_cell_portal"
    }
}
