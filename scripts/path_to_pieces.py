#!/usr/bin/env python


__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2016"
__credits__ = [ "Timothy Tickle", "Brian Haas", "Nir Yosef" ]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@broadinstitute.org"
__status__ = "Development"

#import inspect
import argparse
import os

arg_raw = argparse.ArgumentParser(description="Breaks file paths into pieces and writes them to files.")
arg_raw.add_argument( "--path",
                      dest="path",
                      required = True,
                      help="File path on which to perform manipulations." )
arg_raw.add_argument( "--base_name",
                      dest="base_name",
                      default = None,
                      help="Path to write the base name to (file name no extention)." )
arg_raw.add_argument( "--file_name",
                      dest="file_name",
                      default = None,
                      help="Path to write the file name to (file name with extension devoid of parent directory information)." )
args_parsed = arg_raw.parse_args()

# Base name
if args_parsed.base_name:
    with open(args_parsed.base_name,"w") as base_name_file:
        base_name_file.write(os.path.splitext(os.path.basename(args_parsed.path))[0]+"\n")

# File name
if args_parsed.file_name:
    with open(args_parsed.file_name,"w") as file_name_file:
        file_name_file.write(os.path.basename(args_parsed.path)+"\n")
