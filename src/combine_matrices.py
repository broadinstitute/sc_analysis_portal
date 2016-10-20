#!/usr/bin/env python

"""
This script consolates different (square) matrix files
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

STDOUT = None

class MatrixConsolidator:
    """
    Merges matrix files by row and column name.
    """

    def __init__(self, logfile=STDOUT, delim="\t", NA_value=0):
        self.data = {}
        self.row_order = []
        self.delim = delim
        self.NA = NA_value

        self.logname = "MatrixConsolidator"
        self.logger = logging.getLogger(self.logname)
        if logfile:
            log_settings = logging.FileHandler(filename=logfile, mode="w")
        log_settings.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        self.logger.addHandler(log_settings)
        self.logger.setLevel(logging.INFO)

    def add_matrix(self, matrix_path):
        """ Parse and add square matrix file.
        :param matrix_path: Matrix file to parse and combine
        :type matrix_path: String
        :returns: Indicator of success
        :rtype: Boolean
        """
        header = []
        with open(matrix_path, "r") as in_file:
            in_reader = csv.reader(in_file, delimiter=self.delim)
            for tokens in in_reader:
                # Assume header is the first line and create if we do not have it.
                if not header:
                    header = tokens
                    for header_token in header:
                        if header_token in self.data:
                            self.logger.error(" ".join(["The column name'",
                                                        header_token,
                                                        "'was repeated and could not be added."]))
                            return(False)
                    continue
                row_name = tokens[0]
                self.row_order.append(row_name)
                for column_index in range(1,len(tokens)):
                    # Check to make sure the column is not alreay in the matrix
                    sample_data = self.data.setdefault(header[column_index],{})
                    # Check for an error when repeated row names are encountered
                    if row_name in sample_data:
                        self.logger.error(" ".join(["The row name'",
                                                    row_name,
                                                    "'was repeated in sample'",
                                                    header[column_index]+"'."]))
                        self.row_order = []
                        self.data = {}
                        return(False)
                    sample_data[row_name]=tokens[column_index]
        self.row_order = list(set(self.row_order))
        return(True)

    def write_to_file(self, output_path):
        """ Write output to file.
        :param output_path: Path to write contents of Consolidator
        :type output_path: String
        :returns: Indicator of success
        :rtype: Boolean
        """
        if not output_path or not self.data:
            return(False)
        with open(output_path, "w") as out:
            column_header = list(self.data.keys())
            output = [self.delim.join([""]+column_header)]
            for row_name in self.row_order:
                write_line = [row_name] + [self.data[col_name].get(row_name, self.NA)
                                           for col_name in column_header]
                output.append(self.delim.join(write_line))
            out.write("\n".join(output))
            return(True)

if __name__ == "__main__":

    # Parse arguments
    prsr_arguments = argparse.ArgumentParser(prog='combine_matrix.py',
                                             description='Consolidate matrix files into one matrix output.',
                                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    prsr_arguments.add_argument("--log",
                                metavar="Log File",
                                dest="log_name",
                                default=None,
                                help="Path to the file to use in logging.")

    prsr_arguments.add_argument("--out",
                                metavar="Out File",
                                dest="output_file",
                                default=None,
                                help="Path to the output to create.")

    prsr_arguments.add_argument(dest="in_matrices",
                                nargs="+",
                                help="Matrix file(s) to combine.")
    args = prsr_arguments.parse_args()

    merger = MatrixConsolidator(os.path.abspath(args.log_name))
    success = True
    for tsv_file in args.in_matrices:
        success = success and merger.add_matrix(os.path.abspath(tsv_file))
    if success:
        success = merger.write_to_file(os.path.abspath(args.output_file))
        if success:
            merger.logger.info("Writing output to file "+args.output_file)
        else:
            merger.logger.error("Failed to write to file "+args.output_file)
    else:
        merger.logger.error("Failed to merge tables.")
