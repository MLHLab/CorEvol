#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 14:55:27 2024

@author: abhishake

This script runs alignment over all of the non-redundant core cluster files in temp
"""

import os
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
import json
import argparse


def path_setup(config_file = 'file_paths.json'):
    with open('file_paths.json', 'r') as f:
        paths = json.load(f)

    temp_dir = paths["temp_dir"]
    in_dir = paths["core_clusters"]
    temp_aligned_name = temp_dir+"/Clusters_aligned"
    paths["aligned_core_cluster"] = temp_aligned_name
    return paths, in_dir, temp_aligned_name

def setup_future_paths(paths, temp_aligned_name):
    # paths = {
    #     "temp_dir" : paths["temp_dir"],
    #     "result_dir" : paths["result_dir"],
    #     "fasta_output": paths["fasta_output"],
    #     "csv_output": paths["csv_output"],
    #     "combined_seqs": paths["combined_seqs"],
    #     "core_clusters": paths["core_clusters"],
    #     "rep_seq_list": paths["rep_seq_list"],
    #     "aligned_core_cluster" : temp_aligned_name
    #     }
    paths["core_clusters"]
    paths["aligned_core_cluster"] = temp_aligned_name
    with open('file_paths.json', 'w') as f:
        json.dump(paths, f)


# Loop over each file in the directory
def main(phylogeny_cutoff):
    paths, in_dir, temp_aligned_name = path_setup('file_paths.json')
    os.makedirs(temp_aligned_name)
    # os.makedirs(out_dir)
    for filename in os.listdir(in_dir):
        if filename.endswith('.fa'):  # Process only fasta files
            in_file = os.path.join(in_dir, filename)
            out_file = os.path.join(temp_aligned_name, "aligned_" + filename)

            # Check the length of each sequence in the file
            min_length = float('inf')  # Initialize min_length to infinity
            for record in SeqIO.parse(in_file, 'fasta'):
                seq_length = len(record.seq)
                if seq_length < min_length:
                    min_length = seq_length

            # Run MAFFT only if all sequences are at least 500 nucleotides long
            if min_length >= phylogeny_cutoff:
                # MAFFT command line
                mafft_cline = MafftCommandline(input=in_file)

                # Execute the alignment
                stdout, stderr = mafft_cline()

                # Write the alignment to a file
                with open(out_file, "w") as handle:
                    handle.write(stdout)
            else:
                print(f"Skipping {filename} due to sequence(s) shorter than {phylogeny_cutoff} nucleotides.")
    
    setup_future_paths(paths, temp_aligned_name)


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description="MAFFT alignment")
#     parser.add_argument("-pC", "--phylogeny_cutoff", type=int, default=500, help="Minimum length of necleotides prior to alignment")
#     args = parser.parse_args()
#     main(args)
