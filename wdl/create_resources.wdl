workflow create_resources {

    File in_fasta
    File in_gtf
    File reference_name

    call ref_path_pieces {
        input:file_path=in_fasta, out_dir=reference_name
    }

    call gtf_path_pieces {
        input:file_path=in_gtf, out_dir=reference_name
    }

    call create_all_resources {
        input:fasta=in_fasta, gtf=in_gtf, in_reference_name, ref_base_name = ref_path_pieces.base_name, ref_file_name = ref_path_pieces.file_name, gtf_file_name = gtf_path_pieces.file_name
    }
}

task ref_path_pieces {
    File file_path
    File out_dir
    command <<<
    path_to_pieces.py --path file_path \
                      --base_name "${out_dir}/ref_base_name.txt" \
                      --file_name "${out_dir}/ref_file_name.txt"
    >>>
    output {
        String base_name = read_string("${out_dir}/ref_base_name.txt")
        String file_name = read_string("${out_dir}/ref_file_name.txt")
    }
}

task gtf_path_pieces {
    File file_path
    File out_dir
    command <<<
    path_to_pieces.py --path file_path \
                      --file_name "${out_dir}/gtf_file_name.txt"
    >>>
    output {
        String file_name = read_string("${out_dir}/gtf_file_name.txt")
    }
}

task create_all_resources {

    File fasta
    File gtf
    File reference_name
    String ref_base_name
    String ref_file_name

    command <<<
    prep_resources.py --reference_fasta fasta \
                      --gtf gtf \
                      --out_dir reference_name \
    >>>

    output {
        File out_fasta = "${reference_name}/${ref_file_name}"
        File out_gtf = "${reference_name}/${gtf_file_name}"
        File out_faidx = "${reference_name}/${ref_file_name}.fai"
        File out_dict = "${reference_name}/${ref_base_name}.dict"
        File out_ref_flat = "${reference_name}/${ref_base_name}.refFlat"
        File out_gene_interval = "${reference_name}/${ref_base_name}.genes.intervals"
        File out_rrna_interval = "${reference_name}/${ref_base_name}.rRNA.intervals"
        File out_intergenic_interval = "${reference_name}/${ref_base_name}.intergenic.intervals"
        File rsem_index = "${reference_name}/${ref_base_name}"
    }        
}
