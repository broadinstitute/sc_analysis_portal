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
import sys
# Add sciedpiper
sys.path.insert(0,os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "..", "src", "SciEDPipeR"]))
sys.path.insert(0,os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "..", "src", "SciEDPipeR", "sciedpiper"]))
import sciedpiper.Command as Command
import sciedpiper.PipelineRunner as PipelineRunner


class PrepareSingleCell( PipelineRunner.PipelineRunner ):
    """
    Generates resources for a pipeline for single cell analysis.
    """

    def func_update_arguments( self, arg_raw ):
        """
        Updates to the arg parser, command line options

        * arg_raw : Arguments ( not yet parsed )
                  : Arguments
        * return  : Updated Arguments
                  : Arguments
        """

        arg_raw.prog = "prepare_resources.py"
        arg_raw.description = "Prepare Resources for the Single Cell Analysis Pipeline."
        # Required
        arg_raw.add_argument( "--reference_fasta",
                              dest="reference",
                              required=True,
                              help="Species Reference fasta." )
        arg_raw.add_argument( "--gtf",
                              dest="gtf",
                              required=True,
                              help="Reference GTF." )
        return( arg_raw )


    def func_make_commands( self, args_parsed, cur_pipeline ):
        """
        Generate pipeline commands for analysis.
        """

        # File names
        sample_name = self.func_base_file( args_parsed.reference )
        project_fasta = self.func_update_file_location( args_parsed.reference,
                                                        args_parsed.str_out_dir )
        project_gtf = self.func_update_file_location( args_parsed.gtf,
                                                      args_parsed.str_out_dir )
        rsem_index_prefix = os.path.join( args_parsed.str_out_dir,
                                          sample_name )
        reduced_gtf = self.func_tag_file(project_gtf, "reduced" )
        ref_dict = self.func_switch_ext( project_fasta, ".dict" )
        ref_flat = self.func_switch_ext( project_fasta, ".refFlat" )
        ref_fai = project_fasta + ".fai"
        ref_rrna = self.func_switch_ext( project_fasta, ".rRNA.intervals" )
        ref_gene_interval = self.func_switch_ext( project_fasta, ".genes.intervals" )
        ref_genic_interval = self.func_switch_ext( project_fasta, ".intergenic.intervals" )

        # Make other files given the dependency tree
        commands = []

        # Copy fasta
        cp_fasta_cmdline = " ".join(["cp", args_parsed.reference,
                                     project_fasta])
        cmd_copy_fasta = Command.Command(str_cur_command=cp_fasta_cmdline,
                                  lstr_cur_dependencies=[args_parsed.reference],
                                  lstr_cur_products=[project_fasta])
        commands.append(cmd_copy_fasta)

        # Copy gtf
        cp_gtf_cmdline = " ".join(["cp", args_parsed.gtf,
                                     project_gtf])
        cmd_copy_gtf = Command.Command(str_cur_command=cp_gtf_cmdline,
                                  lstr_cur_dependencies=[args_parsed.gtf],
                                  lstr_cur_products=[project_gtf])
        commands.append(cmd_copy_gtf)

        # Make fasta index
        dict_cmdline = "samtools faidx " + project_fasta
        cmd_fai = Command.Command(str_cur_command=dict_cmdline,
                                  lstr_cur_dependencies=[project_fasta],
                                  lstr_cur_products=[ref_fai])
        commands.append(cmd_fai)

        # Make dict
        dict_cmdline = " ".join(["java -jar picard.jar CreateSequenceDictionary",
                                 "R="+project_fasta,
                                 "O="+ref_dict])
        cmd_dict = Command.Command(str_cur_command=dict_cmdline,
                                     lstr_cur_dependencies=[project_fasta],
                                     lstr_cur_products=[ref_dict])
        commands.append(cmd_dict)

        # Refflat
        refflat_cmdline = " ".join(["ConvertToRefFlat",
                                    "ANNOTATIONS_FILE="+project_gtf,
                                    "SEQUENCE_DICTIONARY="+ref_dict,
                                    "O="+ref_flat])
        cmd_reflat = Command.Command(str_cur_command=refflat_cmdline,
                                     lstr_cur_dependencies=[project_gtf,
                                                            ref_dict],
                                     lstr_cur_products=[ref_flat])
        commands.append(cmd_reflat)

        # GTF
        # Make reduced gtf
        cmd_reduced_gtf = Command.Command(str_cur_command=" ".join(["ReduceGTF",
                                                                    "GTF="+project_gtf,
                                                                    "SEQUENCE_DICTIONARY="+ref_dict,
                                                                    "O="+reduced_gtf]),
                                          lstr_cur_dependencies=[project_gtf, ref_dict],
                                          lstr_cur_products=[reduced_gtf])
        commands.append(cmd_reduced_gtf)

        # Interval files
        cmdline_interval=" ".join(["CreateIntervalsFiles",
                                   "SEQUENCE_DICTIONARY="+ref_dict,
                                   "REDUCED_GTF="+reduced_gtf,
                                   "OUTPUT="+args_parsed.str_out_dir,
                                   "PREFIX="+sample_name])
        cmd_gene_intervals = Command.Command(str_cur_command=cmdline_interval,
                                             lstr_cur_dependencies=[ref_dict,
                                                                    reduced_gtf],
                                             lstr_cur_products=[ref_gene_interval,
                                                                ref_rrna,
                                                                ref_genic_interval])
        commands.append(cmd_gene_intervals)

        ## RSEM Index
        rsem_index_cmdline = " ".join([ "rsem-prepare-reference",
                                        "--bowtie",
                                        "--bowtie2",
                                        "--gtf", project_gtf,
                                        project_fasta,
                                        rsem_index_prefix ])
        rsem_index_deps = [ project_gtf,
                            project_fasta ]
        rsem_index_prods = [ rsem_index_prefix + ".1.bt2",
                             rsem_index_prefix + ".2.bt2",
                             rsem_index_prefix + ".3.bt2",
                             rsem_index_prefix + ".4.bt2",
                             rsem_index_prefix + ".1.ebwt",
                             rsem_index_prefix + ".2.ebwt",
                             rsem_index_prefix + ".3.ebwt",
                             rsem_index_prefix + ".4.ebwt" ]
        rsem_index_cmd = Command.Command( str_cur_command = rsem_index_cmdline,
                                       lstr_cur_dependencies = rsem_index_deps,
                                       lstr_cur_products = rsem_index_prods )
        commands.append( rsem_index_cmd )
        return( commands )


if __name__ == "__main__":

    # Needed to run, calls the script
    PrepareSingleCell().func_run_pipeline()
