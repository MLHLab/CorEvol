#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 17:09:24 2024

@author: abhishake
"""

import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import tempfile
import os
import json
import glob

def get_args():
    """"
    Get the directory containing the FASTA files.

    Returns
    -------
    None.

    """
    desc = "Input the directory containing the FASTA files"
    
    try:
        parser = argparse.ArgumentParser(description=desc)
        parser.add_argument("-i", "--directory", required=True, help="Directory containing the FASTA files")
    
    except NameError:
        sys.stderr.write(
            "An exception occurred with argument parsing. Check your provided options."
        )
        sys.exit(1)
    
    return parser.parse_args()

def main(directory, result):
    # args = get_args()
    output_seqs = []
    original_seqs = []
    headerList = {}
    pre = 1
    column_names = ['original', 'coded']
    df_file_list = pd.DataFrame(columns=column_names)

    # Create a temporary directory
    temp_dir = tempfile.mkdtemp()
    fasta_output_path = os.path.join(temp_dir, 'output.fa')
    csv_output_path = os.path.join(temp_dir, 'header.csv')
    combined_seqs_path = os.path.join(temp_dir, 'combined_seqs.fa')
    filename_map = os.path.join(temp_dir, 'filename_map.csv')
    
    fasta_files = glob.glob(os.path.join(directory, "*"))
    
    os.makedirs(result)

    for fasta_file in fasta_files:
        po = 1
        row_data = {'original': os.path.basename(fasta_file), 'coded': f'{pre:03d}'}
        records = list(SeqIO.parse(fasta_file, "fasta"))
        for record in records:
            new_record = SeqRecord(record.seq)
            new_record.id = f'{pre:03d}' + "_" + f'{po:05d}'
            new_record.description = new_record.id
            headerList[record.id] = new_record.id
            po += 1
            output_seqs.append(new_record)
            original_seqs.append(record)
        row_df = pd.DataFrame([row_data])
        df_file_list = pd.concat([df_file_list, row_df], ignore_index=True)
        pre += 1

    # Write outputs to the temporary directory
    SeqIO.write(output_seqs, fasta_output_path, 'fasta')
    SeqIO.write(original_seqs, combined_seqs_path, 'fasta')
    df = pd.Series(headerList).to_frame()
    df.to_csv(csv_output_path, index=True, header=False)
    df_file_list.to_csv(filename_map, sep = "\t", index=False, header=True)

    # Optionally print out the paths to the output files
    print(f"Temporary directory created: {temp_dir}")
    print(f"FASTA output written to: {fasta_output_path}")
    print(f"Combined original sequences written to: {combined_seqs_path}")
    print(f"CSV output written to: {csv_output_path}")
    paths = {
        "temp_dir" : temp_dir,
        "result_dir" : result,
        "fasta_output": fasta_output_path,
        "csv_output": csv_output_path,
        "combined_seqs": combined_seqs_path,
        "filename_map":filename_map
        }
    
    
    with open('file_paths.json', 'w') as f:
        json.dump(paths, f)


# if __name__ == "__main__":
#     main()