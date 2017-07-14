workflow single_cell_analysis_portal {

    Array[File] sample_names
    File in_biological_batches
    File in_gtf
    File in_index
    File in_left_sample
    File in_negative_controls
    File in_reference
    File in_ref_flat
    File in_ribosome_interval
    File in_right_sample
    File in_technical_batches
    String in_library
    String in_machine
    String in_platform
    String in_sample_name
    Int in_threads

    scatter (name in sample_names){
        call align_qc {
            input: gtf=${in_gtf},
                   reference_fasta=${in_reference},
                   reference_flat=${in_ref_flat},
                   rsem_index=${in_index},
                   rrna_intervals=${in_ribosome_interval},
                   sample_left=${in_left_sample},
                   sample_right=${in_right_sample},
                   library=${in_library},
                   machine=${in_machine},
                   sample_name=${in_sample_name},
                   platform=${in_platform},
                   threads=${in_threads},
        }
        call consolidate_metrics { input: files=run_pipeline.consolidated_qc }
        call consolidate_counts { input: files=run_pipeline.rsem_gene_results }
        call consolidate_tpms { input: files=run_pipeline.rsem_gene_results }
        call scone { input: count_matrix=consolidate_counts.out, technical_batches=${in_technical_batches}, biological_batches=${in_biological_batches}, negative_controls=${in_negative_controls}, qc_matrix=consoliate_metrics.out}
    }

task align_qc {
    File gtf
    File reference_fasta
    File reference_flat
    File rsem_index
    File rrna_intervals
    File sample_left
    File sample_right
    String library
    String machine
    String sample_name
    String platform
    Int threads

    command <<<
        single_cell_analysis.py \
            --reference_fasta ${reference_fasta} \
            --reference_flat ${reference_flat} \
            --rsem_index ${rsem_index} \
            --rrna_intervals ${rrna_intervals} \
            --gtf ${gtf} \
            --sample_left ${sample_left} \
            --sample_right ${sample_right} \
            --add_machine ${machine} \
            --add_name ${sample_name} \
            --add_library ${library} \
            --add_platform ${platform} \
            --threads ${threads} \
    >>>

    output {
        File trim_1_paired = "qc/rnaseqc/${sample_name}_1P.fq.gz"
        File trim_1_unpiared = "qc/rnaseqc/${sample_name}_1U.fq.gz"
        File trim_2_paired = "qc/rnaseqc/${sample_name}_2P.fq.gz"
        File trim_2_unpiared = "qc/rnaseqc/${sample_name}_2U.fq.gz"
        File fastqc = "qc/fast_qc/${sample_name}_fastqc_data.txt"
        File rnaseqc = "qc/rnaseq_qc/${sample_name}_metrics.tsv"
        File rsem_gene_bam = "bam/${sample_name}.genome.bam"
        File rsem_gene_results = "bam/${sample_name}.genes.results"
        File rsem_trans_bam = "bam/${sample_name}.transcript.bam"
        File rsem_trans_results = "bam/${sample_name}.isoforms.results"
        File rsem_stats = "bam/${sample_name}.stats"
        File prepped_bam = "bam/${sample_name}_arg_reorder_dedup.genome.bam"
        File duplicate_metrics = "qc/${sample_name}_duplicate_metrics.txt"
        File rna_seq_metrics = "qc/${sample_name}_collect_rna_seq_metrics.txt"
        File library_metrics = "qc/${sample_name}_estimate_library_complexity.txt"
        File alignment_metrics = "qc/${sample_name}_alignment_summary_metrics.txt"
        File qc_insert_size_txt = "qc/${sample_name}_insert_size.txt"
        File qc_insert_size_pdf = "qc/${sample_name}_insert_size.pdf"
        File consolidated_qc = "qc/${sample_name}_consolidate.tsv"
        File consolidated_qc_log = "qc/${sample_name}_consolidate.log"
    }
}

task consolidate_metrics {
    Array[File] files
    command <<<
        combine_matrix.py \
            --log "qc/qc_matrix.log" \
            --out "qc/qc_matrix.tsv" \
            ${sep=" " files}
    >>>

    output {
        File log = "qc/qc_matrix.log" \
        File out = "qc/qc_matrix.tsv"
    }
}

task consolidate_counts {
    Array[File] files
    File out_dir
    command <<<
        combine_matrix.py \
            --log "raw_counts.tsv" \
            --out "raw_counts.log" \
            --column "Counts" \
            ${sep=" " files}
    >>>
    output {
        File log = "raw_counts.log"
        File out = "raw_counts.tsv"
    }
}

task consolidate_tpms {
    Array[File] files
    command <<<
        combine_matrix.py \
            --log "raw_counts.tsv" \
            --out "raw_counts.log" \
            --column "TPM" \
            ${sep=" " files}
    >>>
    output {
        File log = "raw_tpms.log"
        File out = "raw_tpms.tsv"
    }
}

task scone {
    File count_matrix
    File technical_batches
    File biological_batches
    File negative_controls
    File qc_matrix

    command <<<
        ezbake.R --expr=${count_matrix} \
                 --batch=${technical_batches} \
                 --bio=${biological_batches} \
                 --negcon=${negative_controls} \
                 --qc=${qc_matrix} \
                 --verbose=0 \
                 --norm_adjust_bio=no \
                 --out_dir="normalizations"
    >>>
    output {
        File "normalizations"
    }
}
