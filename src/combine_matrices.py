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

    def __init__(self, logfile=STDOUT, delim="\t", NA_value='0'):
        # Holds all the data
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

    def add_matrix(self, matrix_path, column=None):
        """ Parse and add square matrix file.
        :param matrix_path: Matrix file to parse and combine
        :type matrix_path: String
        :returns: Indicator of success
        :rtype: Boolean
        """
        self.logger.info("Adding matrix=" + matrix_path)
        header = []
        with open(matrix_path, "r") as in_file:
            in_reader = csv.reader(in_file, delimiter=self.delim)
            file_based = None
            target_column_index = None
            if column:
                self.logger.info("Targeting column=" + str(column))
                file_based = "_".join([os.path.splitext(os.path.basename(matrix_path))[0],
                                      column])
            for tokens in in_reader:
                # Assume header is the first line and create if we do not have it.
                if not header:
                    header = tokens
                    # When column is present, reduce header to just that column.
                    if column:
                        if column in tokens:
                            target_column_index = tokens.index(column)
                        else:
                            self.logger.error(" ".join(["The column was not",
                                "found in the matrix.",
                                "Column=",column,
                                "Matrix=",matrix_path]))
                            return(False)
                        if target_column_index == -1:
                            if header_token in self.data:
                                self.logger.error(" ".join(["The following",
                                    "column name could not be found in the",
                                    "matrix. column name =", str(column)]))
                            return(False)
                        header = [header[0], file_based]
                    # Check to make sure the columns have
                    # not already been added in other files.
                    for header_token in header:
                        if header_token in self.data:
                            self.logger.error(" ".join(["The column name'",
                                                        header_token,
                                                        "'was repeated and could not be added."]))
                            return(False)
                    continue
                    self.logger.info("Checked columns and none have already been added and so all are safe to add.")
                # Assumes the first column of the matrix is feature names
                row_name = tokens[0]
                self.row_order.append(row_name)
                # If just using a column reduce to the column
                if column:
                    tokens = [row_name,tokens[target_column_index]]
                for column_index in range(1,len(tokens)):
                    # Check to make sure the column is not already in the matrix
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
            self.logger.error(" ".join(["Could not write to output file.",
                                        "Output_path=", str(output_path),
                                        "Data = ", str(self.data)]))
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

    prsr_arguments.add_argument("--column",
                                metavar="Column_of_interest",
                                dest="target_column",
                                default=None,
                                help="If supplied, only this column will be taken from matrices as they are being combined.")

    prsr_arguments.add_argument("--delim",
                                metavar="Column_of_interest",
                                dest="delim",
                                default="\t",
                                help="Matrix file delimiter.")

    prsr_arguments.add_argument(dest="in_matrices",
                                nargs="+",
                                help="Matrix file(s) to combine.")
    args = prsr_arguments.parse_args()

    if not args.log_name:
        args.log_name = "_".join(["combine"] + [os.path.splitext(os.path.basename(combine_matrix))[0] for combine_matrix in args.in_matrices]) + ".log"
    merger = MatrixConsolidator(logfile=os.path.abspath(args.log_name),
                                delim=args.delim)
    success = True
    for tsv_file in args.in_matrices:
        success = success and merger.add_matrix(matrix_path=os.path.abspath(tsv_file),
                                                column=args.target_column)
    if success:
        success = merger.write_to_file(os.path.abspath(args.output_file))
        if success:
            merger.logger.info("Writing output to file "+args.output_file)
        else:
            merger.logger.error("Failed to write to file "+args.output_file)
    else:
        merger.logger.error("Failed to merge tables.")
