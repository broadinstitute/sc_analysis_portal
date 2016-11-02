#!/usr/bin/env python

"""
This script consolates different types of QC output into a tabular output.
"""

# Import statements
import argparse
import csv
import logging
import os

__author__ = 'Timothy Tickle'
__copyright__ = 'Copyright 2016'
__credits__ = ["Timothy Tickle"]
__license__ = 'BSD-3'
__maintainer__ = 'Timothy Tickle'
__email__ = 'ttickle@broadinstitute.org'
__status__ = 'Development'

class QCParser:
    """
    Parses specific qc files and outputs as a tabular file.
    """

    def __init__(self, logfile):
        self.logname = "QCParserLog"
        self.logger = logging.getLogger(self.logname)
        log_settings = logging.FileHandler(filename=logfile, mode="w")
        log_settings.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        self.logger.addHandler(log_settings)
        self.logger.setLevel(logging.INFO)

    def add_dicts(self, dict_one, dict_two):
        """
        Combines two dicts together.
        Will return false if a key exists in both dicts (error on overwritting).
        Returns a new object, not a reference.

        :param dict_one: Dict to combine
        :type dict_one: Dictionary
        :param dict_two: Dict to combine
        :type dict_two: Dictionary
        :returns: A dict made of the other two merged.
        :rtype: Dict {metric:value}
        """
        dict_return = dict_two.copy()
        for key in dict_one:
            if key in dict_return:
                return(False)
            else:
                dict_return[key]=dict_one[key]
        return(dict_return)

    def consolidate_metrics(self,
                            write_file,
                            sample_name="unidentified_sample",
                            collect_asm_path=None,
                            collect_ism_path=None,
                            collect_rnaseq_path=None,
                            fastq_qc_path=None,
                            estimate_complexity_path=None,
                            rnaseq_qc_path=None):

        collected_metrics = {}

        parse_runs = {estimate_complexity_path : self.parse_estimate_library_complexity,
                      collect_asm_path : self.parse_collect_asm,
                      collect_ism_path : self.parse_collect_ism,
                      collect_rnaseq_path : self.parse_collect_rna_seq_metrics,
                      fastq_qc_path : self.parse_fast_qc,
                      rnaseq_qc_path : self.parse_rna_seq_qc}

        for run_file, run_function in parse_runs.items():
            if run_file:
                run_file = os.path.abspath(run_file)
                if not os.path.exists(run_file):
                    self.logger.error("ERROR, the following file does not exist. File="+run_file)
                    exit(100)
                info = run_function(run_file)
                collected_metrics = self.add_dicts(info, collected_metrics)
                if not collected_metrics:
                    self.logger.error("ERROR, could not parse file. File="+run_file)
                    exit(101)

        self.write_info_to_file(data=collected_metrics,
                                sample_name=sample_name,
                                output_file=write_file)

    def make_safe_key(self, key):
        return("".join([letter if letter.isalnum() else "_" for letter in key]))

    def parse_fast_qc(self, fastqc_file):
        """ Parse FastQC output file.
        :param fastqc_file: File to parse
        :type fastqc_file: String
        :returns: A dict of sample to metric pairings (sample name = key)
        :rtype: Dict {metric:value}
        """
        #TODO adapter contamination.
        # Stored for output.
        output ={}
        dup_key = "#Total Deduplicated Percentage"
        with open(fastqc_file, "r") as qc_file:
            qc_reader = csv.reader(qc_file)
            for tokens in qc_reader:
                if tokens[0] == dup_key:
                    output[self.make_safe_key(dup_key)] = tokens[1]
                    return(output)
        return(output)

    def parse_rna_seq_qc(self, RNASeqQC_file):
        """ Parse RNA-SeQC.jar output file.
        :param RNASeqQC_file: File to parse
        :type RNASeqQC_file: String
        :returns: A dict of sample to metric pairings (sample name = key)
        :rtype: Dict {metric:value}
        """
        #TODO what do we need from there
        # (1) number of reads
        # (2) number of aligned reads
        # (3) percentage of aligned reads
        # (4) percentage of transcripts identified (compared with the
        ## overall number of transcripts identified by at least
        ## one cell in a given run),
        # (5) percentage of duplicate reads
        # (6) primer sequence contamination
        # (7) insert size (mean)
        # (8) insert size (std)
        # (9) complexity
        # (10) percentage of ribosomal reads
        # (11) percentage of coding reads
        # (12) percentage of UTR reads
        # (13) percentage of Intronic reads
        # (14) percentage of intergenic reads
        # (15) percentage of mRNA reads
        # (16) the coefficient of variation of coverage
        # (17) mean 5 Bias
        # (18) mean 3 Bias
        # (19) mean 5 to 3 Bias.
        return(self.parse_piccard_output(in_file=RNASeqQC_file,
                                         header_key="Sample",
                                         metrics=[]))

    def parse_collect_rna_seq_metrics(self, CollectRNAMetrics_file):
        """ Parse CollectRnaSeqMetrics.jar output file.
        :param CollectRNAMetrics_file: File to parse
        :type CollectRNAMetrics_file: String
        :returns: A dict of sample to metric pairings (sample name = key)
        :rtype: Dict {metric:value}
        """
        # TODO
        return({})

    def parse_estimate_library_complexity(self, estimateLC_file):
        """ Parse EstimateLibraryComplexity.jar output.
        :param EstimateLC_file: File to parse
        :type EstimateLC_file: String
        :returns: A dict of sample to metric pairings (sample name = key)
        :rtype: Dict {metric:value}
        """
        # TODO is READ_PAIR_DUPLICATES correct?
        return(self.parse_piccard_output(in_file=estimateLC_file,
                                         header_key="LIBRARY",
                                         metrics=["READ_PAIR_DUPLICATES",
                                                  "ESTIMATED_LIBRARY_SIZE"]))

    def parse_piccard_output(self, in_file, header_key, metrics):
        # Stored for output.
        output ={}
        # Header and index
        header_metrics = None
        with open(in_file, "r") as qc_file:
            qc_reader = csv.reader(qc_file,delimiter="\t")
            for tokens in qc_reader:
                if not tokens or tokens[0][0]=="#":
                    continue
                if tokens[0] == header_key and not header_metrics:
                    header_metrics = dict(zip(tokens,range(len(tokens))))
                    continue
                # Check metric.
                for save_metric in metrics:
                    output[self.make_safe_key(save_metric)] = tokens[header_metrics[save_metric]]
                return(output)
        return(output)

    def parse_collect_asm(self, collectASM_file):
        """ Parse CollectAlignmentSummaryMetrics.jar output.
        :param collectASM_file: File to parse
        :type collectASM_file: String
        :returns: A dict of sample to metric pairings (sample name = key)
        :rtype: Dict {metric:value}
        """
        metric_class = ["PAIR"]
        # Metrics: "NREADS, RALIGN, and NALIGN"
        metrics = ["TOTAL_READS","PF_READS_ALIGNED"]
        # Stored for output.
        output ={}
        # Header and index
        header_key = "CATEGORY"
        header_metrics = None
        with open(collectASM_file, "r") as qc_file:
            qc_reader = csv.reader(qc_file, delimiter="\t")
            for tokens in qc_reader:
                if not tokens or tokens[0][0]=="#":
                    continue
                # Check metrics class.
                if tokens[0] in metric_class:
                    # Check metric.
                    for save_metric in metrics:
                        output[self.make_safe_key(save_metric)] = tokens[header_metrics[save_metric]]
                    continue
                if tokens[0] == header_key and not header_metrics:
                    header_metrics = dict(zip(tokens,range(len(tokens))))
                    continue
        return(output)

    def parse_collect_ism(self, collectISM_file):
        """ Parse CollectInsertSizeMetrics.jar output.
        :param CollectISM_file: File to parse
        :type CollectISM_file: String
        :returns: A dict of sample to metric pairings (sample name = key)
        :rtype: Dict {metric:value}
        """
        return(self.parse_piccard_output(in_file=collectISM_file,
                                      header_key="MEDIAN_INSERT_SIZE",
                                      metrics=["MEAN_INSERT_SIZE",
                                               "STANDARD_DEVIATION"]))

    def write_info_to_file(self, data, sample_name, output_file):
        self.logger.info("Output written to file="+output_file)
        with open(output_file, "w") as out:
            out.write("".join(["metric"+"\t"+sample_name+"\n"]+[key+"\t"+value+"\n" for key, value in data.items()]))

if __name__ == "__main__":

    # Parse arguments
    prsr_arguments = argparse.ArgumentParser(prog='consolidate_qc.py',
                                             description='Consolidate QC Reports in to a Tabular Output.',
                                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    prsr_arguments.add_argument("--out_file",
                                metavar="Output File",
                                dest="output_file",
                                required=True,
                                help="Output File to write to.")
    prsr_arguments.add_argument("--log_file",
                                metavar="Log File",
                                dest="log_file",
                                required=True,
                                help="Log File to write to.")
    prsr_arguments.add_argument("--sample_name",
                                metavar="Sample Name",
                                dest="sample",
                                required=True,
                                help="The name of the sample.")
    prsr_arguments.add_argument("--collect_alignment",
                                metavar="Collect_Alignment_Summary_Metrics_File",
                                dest="collect_asm_path",
                                default=None,
                                help="Path to the output of CollectAlignmentSummaryMetrics.jar.")
    prsr_arguments.add_argument("--collect_insert",
                                metavar="Collect_Insert_Size_Metrics_File",
                                dest="collect_ism_path",
                                default=None,
                                help="Path to the output of CollectInsertSizeMetrics.jar.")
    prsr_arguments.add_argument("--collect_rnaseq",
                                metavar="Collect_RnaSeq_Metrics_File",
                                dest="collect_rnaseq_path",
                                default=None,
                                help="Path to the output of CollectRnaSeqMetrics.jar.")
    prsr_arguments.add_argument("--fastqqc",
                                metavar="FastqQC_File",
                                dest="fastq_qc_path",
                                default=None,
                                help="Path to the output of FastqQC.jar.")
    prsr_arguments.add_argument("--complexity",
                                metavar="Estimate_Library_Complexity_File",
                                dest="estimate_complexity_path",
                                default=None,
                                help="Path to the output of EstimateLibraryComplexity.jar.")
    prsr_arguments.add_argument("--rnaseqqc",
                                metavar="RNASeqQC_File",
                                dest="rnaseq_qc_path",
                                default=None,
                                help="Path to the output of RNASeqQC.jar.")

    args = prsr_arguments.parse_args()

    QCParser(args.log_file).consolidate_metrics(write_file=args.output_file,
                                   sample_name=args.sample,
                                   collect_asm_path=args.collect_asm_path,
                                   collect_ism_path=args.collect_ism_path,
                                   collect_rnaseq_path=args.collect_rnaseq_path,
                                   fastq_qc_path=args.fastq_qc_path,
                                   estimate_complexity_path=args.estimate_complexity_path,
                                   rnaseq_qc_path=args.rnaseq_qc_path)
