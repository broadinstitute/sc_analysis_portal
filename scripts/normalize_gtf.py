#!/usr/bin/env python

"""
Insert gene name if it doesn't exist for a feature.
"""


import argparse
import os
import sciedpiper.Command as Command
import sciedpiper.ParentScript as ParentScript
import csv

__author__ = [ "Asma Bankapur", "Timothy Tickle" ]
__copyright__ = "Copyright 2016"
__credits__ = ["Asma Bankapur", "Timothy Tickle", "Brian Haas"]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@broadinstitute.org"
__status__ = "Development"

GENE_NAME = "gene_name"
GENE_ID = "gene_id"
NAME = "name"
TRANS_NAME = "transcript_name"
TRANS_ID = "transcript_id"
UNKNOWN = '"unknown"'
GTF_COMMENT = "#!"
DELIMITER = "\t"
SPLIT_TOKENS = "; "
SPLIT_COUPLE = " "
STRIP_QUOTES = '"'
DEF_GEN_NAME = ""
DEF_GEN_ID = ""
JOIN_NAID = "_"
DEF_TRANS_NAME = ""
DEF_TRANS_ID = ""
TOKEN_SEP = ";"

def insert_gene_name( gtf_file ):
    corrected_line = [ ]
    with open( gtf_file , "r" ) as gtf:
        gtf_reader = csv.reader(gtf ,delimiter=DELIMITER)
        for feature_row in gtf_reader:
            ##Add comment feature as is and continue
            if feature_row[ 0 ].startswith( GTF_COMMENT ):
                corrected_line.append(DELIMITER.join(feature_row))
                continue

            ###
            #genid: dict for last element in feature_row
            #gene_id_order_list: to maintain last element token order
            ###
            genid = {}
            gene_id_order_list = []

            gene_cood_chr_info = feature_row[ 0:-1 ]
            token_list_coupled = feature_row[-1].split(SPLIT_TOKENS)

            #Loop to split on each token pair and append to order list
            #Also store the token pair in genid as key val
            for token_couple in token_list_coupled:
                token, token_value = token_couple.split(SPLIT_COUPLE)
                gene_id_order_list.append(token)
                genid[token] = token_value

            #Get 'gene_name', if not empty str and strip ""
            #If not new_name then get 'name' and strip ""
            #Get gene_id else empty str and strip
            #If gene_id is not empty then concat name and id
            #If new_name is _ then place UNKNOWN else keep new_name with quotes
            #add 'gene_name' to dict, join with '; ' token couple
            #reconstruct with \t
            new_name = genid.get(GENE_NAME, DEF_GEN_NAME).strip(STRIP_QUOTES)
            new_name = new_name if new_name else genid.get(NAME, DEF_GEN_NAME).strip(STRIP_QUOTES)
            gene_id = genid.get(GENE_ID, DEF_GEN_ID).strip(STRIP_QUOTES)
            new_name = JOIN_NAID.join([id_token for id_token in [new_name, gene_id] if id_token])
            new_name = UNKNOWN if not new_name else STRIP_QUOTES+new_name+STRIP_QUOTES

            trans_id = genid.get(TRANS_ID, DEF_TRANS_ID).strip(STRIP_QUOTES)
            new_transname = genid.get(TRANS_NAME, DEF_TRANS_NAME).strip(TOKEN_SEP).strip(STRIP_QUOTES)
            new_transname = JOIN_NAID.join([id_token for id_token in [new_transname, trans_id] if id_token])
            new_transname = UNKNOWN if not new_transname else STRIP_QUOTES+new_transname+STRIP_QUOTES
            genid[GENE_NAME] = new_name
            genid[TRANS_NAME] = new_transname
            for id_token in [GENE_NAME,TRANS_NAME]:
                if id_token not in gene_id_order_list:
                    gene_id_order_list.append(id_token)
            gene_id_str = SPLIT_TOKENS.join([token+SPLIT_COUPLE+(genid[token]).strip(TOKEN_SEP)
                                     for token in gene_id_order_list])
            gene_id_str += TOKEN_SEP
            corrected_line.append(DELIMITER.join(feature_row[0:-1] + [gene_id_str]))
    return(corrected_line)

def write_to_file( recontructed_feature_list, modified_gtf_file ):
    with open( modified_gtf_file , "w" ) as mgh:
         mgh.write( "\n".join(recontructed_feature_list))

if __name__ == "__main__":
   arg_raw = argparse.ArgumentParser(prog="normalize_gtf.py",description="Normalize GTF to well-mannered formats.")
   arg_raw.add_argument( "--gtf_file", help="GTF file for normalizing." )
   arg_raw.add_argument( "--output_gtf_path", help="Output GTF file path." )
   args_parsed = arg_raw.parse_args( )
   recontructed_feature_list = insert_gene_name( args_parsed.gtf_file )
   write_to_file( recontructed_feature_list, args_parsed.output_gtf_path )
