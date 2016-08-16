#!/usr/bin/env python


__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2016"
__credits__ = [ "Timothy Tickle", "Brian Haas", "Nir Yosef" ]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@broadinstitute.org"
__status__ = "Development"

#import inspect
import os
import sciedpiper.Command as Command
import sciedpiper.PipelineRunner as PipelineRunner

# Directories
C_BAM = "bam"
C_QC = "qc"

class SingleCell( PipelineRunner.PipelineRunner ):
    """
    Generates a pipeline for single cell analysis.
    This includes a wide swath of analysis which can
    be turned on or off as needed.
    """

    def func_update_arguments(self, arg_raw ):
        """
        Updates to the arg parser, command line options

        * arg_raw : Arguments ( not yet parsed )
                  : Arguments
        * return  : Updated Arguments
                  : Arguments
        """

        arg_raw.prog = "single_cell_analysis.py"
        arg_raw.description = "Single Cell Analysis Pipeline."
        # Required
        arg_raw.add_argument("--no_trim",
                             dest="do_trim",
                             default=True,
                             action="store_false",
                             help="When supplied, turns off trimmmomatic." )
        arg_raw.add_argument("--reference_fasta",
                             dest="reference",
                             required=True,
                             help="Species Reference fasta." )
        arg_raw.add_argument("--reference_flat",
                             dest="reference_flat",
                             required=True,
                             help="Flat file for species reference." )
        arg_raw.add_argument("--rsem_index",
                             dest="rsem_index",
                             required=True,
                             help="RSEM index for reference." )
        arg_raw.add_argument("--rrna_intervals",
                             dest="rrna_intervals",
                             required=True,
                             help="Interval file for ribosomal sequences." )
        arg_raw.add_argument("--gtf",
                             dest="gtf",
                             required=True,
                             help="Reference GTF." )
        arg_raw.add_argument("--sample_left",
                             dest="sample_left",
                             required=True,
                             help="Left Sample fastq.gz." )
        arg_raw.add_argument("--sample_right",
                             dest="sample_right",
                             required=True,
                             help="Right Sample fastq.gz." )
        # Optional
        arg_raw.add_argument("--add_machine",
                             dest="machine",
                             default="machine",
                             help="Machine name to be added to bam." )
        arg_raw.add_argument("--add_name",
                             dest="name",
                             default="name",
                             help="Sample name to be added to bam." )
        arg_raw.add_argument("--add_library",
                             dest="library",
                             default="library",
                             help="Library name to be added to bam." )
        arg_raw.add_argument("--add_plateform",
                             dest="plateform",
                             default="plateform",
                             help="Plateform name to be added to bam." )
        arg_raw.add_argument("--threads",
                             dest="threads",
                             default=1,
                             help="Number of threads for underlying tools." )
        return( arg_raw )


    def func_make_commands( self, args_parsed, cur_pipeline ):
        """
        Generate pipeline commands for analysis.
        """

        # File names
        alignment_summary_file = os.path.join( args_parsed.str_out_dir,
                                               C_QC,
                                               "alignment_summary_metrics.txt" )
        collect_rna_metric_file = os.path.join( args_parsed.str_out_dir,
                                                C_QC,
                                                "collect_rna_seq_metrics.txt" )
        complexity_file = os.path.join( args_parsed.str_out_dir,
                                        C_QC,
                                        "estimate_library_complexity.txt" )
        fastqc_out_dir = os.path.join( args_parsed.str_out_dir,
                                       C_QC,
                                       "fastqc" )
        in_fastq_left = args_parsed.sample_left
        in_fastq_right = args_parsed.sample_right
        in_trimmed_left = self.func_tag_file( in_fastq_left ,"trim" )
        in_trimmed_right = self.func_tag_file( in_fastq_right,"trim" )
        insert_size_summary_file = os.path.join( args_parsed.str_out_dir,
                                                 C_QC,
                                                 "insert_size.txt" )
        insert_size_summary_pdf = os.path.join( args_parsed.str_out_dir,
                                                C_QC,
                                                "insert_size.pdf" )
        mark_metrics_file = os.path.join( args_parsed.str_out_dir,
                                          C_QC,
                                          "duplicates_metrics.txt" )
        rnaseqc_out_dir = os.path.join( args_parsed.str_out_dir,
                                        C_QC,
                                        "rnaseqc" )
        rsem_output = os.path.join( args_parsed.str_out_dir,
                                    C_BAM )
        sample_name = self.func_base_file( args_parsed.str_out_dir )
        # bams
        add_rg_bam = self.func_tag_file( rsem_output, "arg" )
        reorder_bam = self.func_tag_file( add_rg_bam, "reorder" )
        dedup_bam = self.func_tag_file( reorder_bam, "dedup" )

        # Make other files given the dependency tree
        commands = []

        # Trimmomatic
        if args_parsed.do_trim:
            paired_forward_out = os.path.join( rnaseqc_out_dir,
                                              sample_name+"_1P.fq.gz" )
            unpaired_forward_out = os.path.join( rnaseqc_out_dir,
                                                sample_name+"_1U.fq.gz" )
            paired_reverse_out = os.path.join( rnaseqc_out_dir,
                                              sample_name+"_2P.fq.gz" )
            unpaired_reverse_out = os.path.join( rnaseqc_out_dir,
                                                sample_name+"_2U.fq.gz" )
            trim_cmdline = " ".join([ "java -jar trimmomatic.jar PE",
                                      "-threads", str(args_parsed.threads),
                                      "-trimlog", rnaseqc_out_dir+"_log.txt",
                                      in_fastq_left,
                                      in_fastq_right,
                                      paired_forward_out,
                                      unpaired_forward_out,
                                      paired_reverse_out,
                                      unpaired_reverse_out,
                                      "LEADING:15",
                                      "TRAILING:15",
                                      "MINLEN:36",
                                      "2>",
                                      os.path.join( rnaseqc_out_dir,
                                                    "trimmomatic.log.stats" )])
            trim_deps = [ in_fastq_left,
                          in_fastq_right ]
            trim_prods = [ paired_forward_out,
                           unpaired_forward_out,
                           paired_reverse_out,
                           unpaired_reverse_out ]
            trim_cmd = Command.Command( str_cur_command = trim_cmdline,
                                        lstr_cur_dependencies = trim_deps,
                                        lstr_cur_products = trim_prods )
            commands.append( trim_cmd )

        ## FastQC
        fastqc_cmdline = " ".join([ "fastqc",
                                    in_fastq_left,
                                    in_fastq_right,
                                    "--outdir="+fastqc_out_dir ])
        fastqc_deps = [ in_fastq_left, in_fastq_right ]
        fastqc_prods = [ fastqc_out_dir ]
        fastqc_cmd = Command.Command( str_cur_command = fastqc_cmdline,
                                      lstr_cur_dependencies = fastqc_deps,
                                      lstr_cur_products = fastqc_prods )
        commands.append( fastqc_cmd )

        # RSEM
        rsem_cmdline = " ".join([ "rsem-calculate-expression",
                                  "-p", str(args_parsed.threads),
                                  "--paired-end",
                                  "--bowties2",
                                  "--estimate-rspd",
                                  "--output_genome-bam",
                                  in_fastq_left,
                                  in_fastq_right,
                                  args_parsed.rsem_index,
                                  rsem_output ])
        rsem_deps = [ in_fastq_left,
                      in_fastq_right,
                      args_parsed.rsem_index ]
        rsem_prods = [ rsem_output ]
        rsem_cmd = Command.Command( str_cur_command = rsem_cmdline,
                                    lstr_cur_dependencies = rsem_deps,
                                    lstr_cur_products = rsem_prods )
        commands.append( rsem_cmd )

        # Prep for QC
        ## Add_rg
        add_rg_cmdline = " ".join([ "java -jar AddOrReplaceReadGroups.jar",
                                    "I="+rsem_output,
                                    "O="+add_rg_bam,
                                    "SO=coordinate",
                                    "RGLB="+args_parsed.library,
                                    "RGPL="+args_parsed.plateform,
                                    "RGPU="+args_parsed.machine,
                                    "RGSM="+args_parsed.name ])
        add_rg_deps = [ rsem_output ]
        add_rg_prods = [ add_rg_bam ]
        add_rg_cmd = Command.Command( str_cur_command = add_rg_cmdline,
                                      lstr_cur_dependencies = add_rg_deps,
                                      lstr_cur_products = add_rg_prods )
        commands.append( add_rg_cmd )
        ## Reorder_bam
        reorder_cmdline = " ".join([ "java -Xmx4g -jar ReorderSam.jar",
                                     "I="+add_rg_bam,
                                     "R="+args_parsed.reference,
                                     "O="+reorder_bam ])
        reorder_deps = [ add_rg_bam, args_parsed.reference ]
        reorder_prods = [ reorder_bam ]
        reorder_cmd = Command.Command( str_cur_command = reorder_cmdline,
                                       lstr_cur_dependencies = reorder_deps,
                                       lstr_cur_products = reorder_prods )
        commands.append( reorder_cmd )
        ## Mark Duplicates
        dup_cmdline = " ".join([ "java -jar MarkDuplicates.jar",
                                 "I="+reorder_bam,
                                 "O="+dedup_bam,
                                 "CREATE_INDEX=true",
                                 "M="+mark_metrics_file ])
        dup_deps = [ reorder_bam ]
        dup_prods = [ dedup_bam, mark_metrics_file ]
        dup_cmd = Command.Command( str_cur_command = dup_cmdline,
                                   lstr_cur_dependencies = dup_deps,
                                   lstr_cur_products = dup_prods )
        commands.append( dup_cmd )

        # Run QC
        ## RNASeqQC
        rnaseqc_cmdline = " ".join([ "java -Xmx4g -jar RNA-SeQC.jar",
                                     "-n 1000",
                                     "-s \"sample|\""+dedup_bam+"\"|RNASeQC",
                                     "-t", args_parsed.gtf,
                                     "-r", args_parsed.reference,
                                     "-o", rnaseqc_out_dir])
        rnaseqc_deps = [ dedup_bam,
                         args_parsed.gtf,
                         args_parsed.reference ]
        rnaseqc_prods = [ rnaseqc_out_dir ]
        rnaseqc_cmd = Command.Command( str_cur_command = rnaseqc_cmdline,
                                       lstr_cur_dependencies = rnaseqc_deps,
                                       lstr_cur_products = rnaseqc_prods )
        commands.append( rnaseqc_cmd )
        ## Collect RNASeq Metrics
        rnametrics_cmdline = " ".join([ "java -jar CollectRnaSeqMetrics.jar",
                                        "I="+dedup_bam,
                                        "O="+collect_rna_metric_file,
                                        "REF_FLAT="+args_parsed.reference_flat,
                                        "STRAND=NONE",
                                        "RIBOSOMAL_INTERVALS="+args_parsed.rrna_intervals ])
        rnametrics_deps = [ dedup_bam,
                            args_parsed.reference_flat,
                            args_parsed.rrna_intervals ]
        rnametrics_prods = [ collect_rna_metric_file ]
        rnametrics_cmd = Command.Command( str_cur_command = rnametrics_cmdline,
                                       lstr_cur_dependencies = rnametrics_deps,
                                       lstr_cur_products = rnametrics_prods )
        commands.append( rnametrics_cmd )
        ## Estimate Library Complexity
        est_cmdline = " ".join([ "java -jar EstimateLibraryComplexity.jar",
                                 "I="+dedup_bam,
                                 "O="+complexity_file ])
        est_deps = [ dedup_bam ]
        est_prods = [ complexity_file ]
        est_cmd = Command.Command( str_cur_command = est_cmdline,
                                   lstr_cur_dependencies = est_deps,
                                   lstr_cur_products = est_prods )
        commands.append( est_cmd )
        ## Collect Alignment Summary Metrics
        algnmetric_cmdline = " ".join([ "java -jar CollectAlignmentSummaryMetrics.jar",
                                        "REFERENCE_SEQUENCE="+args_parsed.reference,
                                        "INPUT="+dedup_bam,
                                        "OUTPUT="+alignment_summary_file ])
        algnmetric_deps = [ dedup_bam,
                            args_parsed.reference ]
        algnmetric_prods = [ alignment_summary_file ]
        algnmetric_cmd = Command.Command( str_cur_command = algnmetric_cmdline,
                                        lstr_cur_dependencies = algnmetric_deps,
                                        lstr_cur_products = algnmetric_prods )
        commands.append( algnmetric_cmd )
        ## Collect Insert Size
        insert_cmdline = " ".join([ "java -jar CollectInsertSizeMetrics.jar",
                                    "I="+dedup_bam,
                                    "O="+insert_size_summary_file,
                                    "H="+insert_size_summary_pdf,
                                    "M=0.5" ])
        insert_deps = [ dedup_bam ]
        insert_prods = [ insert_size_summary_file, insert_size_summary_pdf ]
        insert_cmd = Command.Command( str_cur_command = insert_cmdline,
                                      lstr_cur_dependencies = insert_deps,
                                      lstr_cur_products = insert_prods )
        commands.append( insert_cmd )

        # Prep and run for Scone (Normalization)
        ## Collect Metrics for Scone
        #prepscone_cmdline = " ".join()
        #prepscone_deps = []
        #prepscone_prods = []
        #prepscone_cmd = Command.Command( str_cur_command = prepscone_cmdline,
        #                                 lstr_cur_dependencies = prepscone_deps,
        #                                 lstr_cur_products = prepscone_prods )
        #commands.append( prepscone_cmd )
        ## Scone
        #scone_cmdline = " ".join()
        #scone_deps = []
        #scone_prods = []
        #scone_cmd = Command.Command( str_cur_command = scone_cmdline,
        #                             lstr_cur_dependencies = scone_deps,
        #                             lstr_cur_products = scone_prods )
        #commands.append( scone_cmd )
        return commands


if __name__ == "__main__":

    # Needed to run, calls the script
    SingleCell().func_run_pipeline()
