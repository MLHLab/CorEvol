#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 14:00:56 2023

@author: abhishake
perform cd-hit, find core clusters, removes paralogs and takes the highest identity entries among them, generates individual fasta files for each core cluster.
"""
from pycdhit import read_fasta, CDHIT
from Bio import SeqIO
import os
import json
import argparse
import pandas as pd

def run_cdhit(input_path, cdhit_path="/usr/bin", c=0.9, d=0, sc=1, aL=0.9, aS=0.9, g=1, n=9):
    """
    Run the CD-HIT clustering algorithm on the input sequence data.

    This function takes an input FASTA file containing sequences and performs clustering using the CD-HIT algorithm.
    It returns two DataFrames: one containing the representative sequences of each cluster and another containing the
    cluster membership information.

    Parameters:
        input_path (str): The path to the input FASTA file containing the sequences to be clustered.
        cdhit_path (str, optional): The path to the CD-HIT executable. Default is "/usr/bin".
        c (float, optional): The sequence identity threshold. Sequences with identity higher than this threshold will be
            clustered together. Default is 0.9 (90% identity).
        d (int, optional): The length difference cutoff. If set to a positive value, sequences shorter than this value
            will be removed. Default is 0 (no length difference cutoff).
        sc (int, optional): The length difference cutoff in amino acid for sequences to be clustered. Default is 1.
        aL (float, optional): The alignment coverage for the longer sequence. Default is 0.9 (90% coverage).
        aS (float, optional): The alignment coverage for the shorter sequence. Default is 0.9 (90% coverage).
        g (int, optional): The maximum available memory in GB. Default is 1 GB.
        n (int, optional): The word length. Default is 3.

    Returns:
        tuple: A tuple containing two DataFrames:
            - df_out (pandas.DataFrame): A DataFrame containing the representative sequences of each cluster.
                The DataFrame has the following columns:
                    - "cluster_id": The unique identifier of the cluster.
                    - "sequence_id": The identifier of the representative sequence.
                    - "sequence": The representative sequence itself.
            - df_clstr (pandas.DataFrame): A DataFrame containing the cluster membership information.
                The DataFrame has the following columns:
                    - "cluster_id": The unique identifier of the cluster.
                    - "sequence_id": The identifier of the sequence belonging to the cluster.

    Example:
        # Run CD-HIT clustering on an input FASTA file
        input_file = "sequences.fasta"
        df_out, df_clstr = run_cdhit(input_file, c=0.95, d=10, aL=0.8, aS=0.8)

        # Print the representative sequences
        print(df_out)

        # Print the cluster membership information
        print(df_clstr)

    Notes:
        - The function requires the `pycdhit` library to be installed.
        - The input FASTA file should contain sequences in the standard FASTA format.
        - The CD-HIT executable should be properly installed and accessible from the specified `cdhit_path`.
        - The function returns two DataFrames: `df_out` containing the representative sequences of each cluster, and
          `df_clstr` containing the cluster membership information.
        - Adjust the clustering parameters (`c`, `d`, `sc`, `aL`, `aS`, `g`, `n`) according to your specific requirements.

    See also:
        - `pycdhit.read_fasta`: Function to read sequences from a FASTA file.
        - `pycdhit.CDHIT`: Class representing the CD-HIT clustering algorithm.
    """
    df_in = read_fasta(input_path)
    cdhit = CDHIT(prog="cd-hit-est", path=cdhit_path)
    df_out,df_clstr = cdhit.set_options(c=c, d=d, sc=sc, aL=aL, aS=aS, g=g, n=n).cluster(df_in)
    return df_out, df_clstr



def path_setup(config_file = 'file_paths.json'):
    with open('file_paths.json', 'r') as f:
        paths = json.load(f)

    temp_dir = paths["temp_dir"]
    rep_seq_path = os.path.join(temp_dir, 'representative_subset.csv')
    master_clstr_path = os.path.join(paths["result_dir"],'Cluster_matrix.csv')
    paths["core_clusters"] = temp_dir+"/Clusters"
    temp_clus_dir = os.makedirs(temp_dir+"/Clusters")
    return paths, rep_seq_path, master_clstr_path, temp_clus_dir


def filter_fasta_file(input_file, identifiers, cluster, temp_clus_name):
    cluster_file_path = os.path.join(temp_clus_name, f'{cluster}.fa')
    with open(input_file, 'r') as fasta_file, open(cluster_file_path, 'w') as output:
        sequences = SeqIO.parse(fasta_file, 'fasta')
        filtered = (sequence for sequence in sequences if sequence.id in identifiers)
        SeqIO.write(filtered, output, 'fasta')
        
def setup_future_paths(paths, rep_seq_path):
    # paths = {
    #     "temp_dir" : paths["temp_dir"],
    #     "result_dir" : paths["result_dir"],
    #     "fasta_output": paths["fasta_output"],
    #     "csv_output": paths["csv_output"],
    #     "combined_seqs": paths["combined_seqs"],
    #     "core_clusters": paths["temp_clus_name"],
    #     "rep_seq_list": rep_seq_path,
    #     "master_file" : master_clstr_path,
    #     "filename_map":paths["filename_map"]
    #     }
    # paths["temp_clus_name"] = core_clusters
    paths["rep_seq_list"] = rep_seq_path
    with open('file_paths.json', 'w') as f:
        json.dump(paths, f)        
        
def map_values(x, value_map):
    if pd.isna(x):
        return x
    else:
        values = str(x).split(',')
        mapped_values = [value_map.get(value.strip(), value.strip()) for value in values]
        return ','.join(mapped_values)

def transform_clstr2mas(df_clstr, header_df, filename_map):
    # Create a new DataFrame using pivot_table
    new_df = df_clstr.pivot_table(index='cluster', columns='Org', values='identifier', aggfunc=lambda x: ', '.join(x))
         
    # Create a dictionary mapping coded headers to actual sequence names from header_df
    header_map = dict(zip(header_df['coded_seq_name'], header_df['actual_seq_name']))
    
    # # Create a dictionary mapping organism codes to actual organism names from filename_map
    filename_map['original'] = filename_map['original'].apply(lambda x: os.path.splitext(x)[0])
    org_map = dict(zip(filename_map['coded'], filename_map['original']))
    
    # # Replace the organism codes in the column names of new_df with actual organism names using the org_map
    new_mas_df = new_df.rename(columns=org_map)
    
    # # Replace the coded entries in new_mas_df with actual sequence names using the header_map
    new_mas_df = new_mas_df.applymap(lambda x: map_values(x, header_map))
    
    # Fill any missing values with an empty string
    new_mas_df.fillna('', inplace=True)
    
    return new_mas_df

def main(path="/usr/bin", identity=0.9, length_diff=0, length_cutoff=1, align_cov_long=0.9, align_cov_short=0.9, memory=4, word_length=9):

    
    paths, rep_seq_path, master_clstr_path, temp_clus_dir = path_setup('file_paths.json') 
           
    df_out, df_clstr = run_cdhit(paths["fasta_output"], cdhit_path= path, c=identity, d=length_diff, sc=length_cutoff,aL=align_cov_long, aS=align_cov_short, g=memory, n=word_length)

    df_rep = df_clstr[df_clstr['is_representative'] == True]
    df_rep = df_clstr[df_clstr['is_representative'] == True].copy()  

    df_rep.loc[:, 'Cluster_ID'] = 'Cluster_' + df_rep['cluster'].astype(str).str.zfill(3)

    df_rep.to_csv(rep_seq_path, sep = "\t", index=False)

    df_clstr[["Org","SeqCount"]] = df_clstr["identifier"].str.split("_",expand=True)
    
    header_df = pd.read_csv(paths['csv_output'], names=["actual_seq_name","coded_seq_name"])
    
    filename_map = pd.read_csv(paths['filename_map'], sep = "\t", dtype=str, header=0)
    mas_df = transform_clstr2mas(df_clstr, header_df,filename_map)
    
    mas_df.to_csv(master_clstr_path,sep = "\t")

    df_grouped = df_clstr.groupby(['cluster', 'Org']).size().unstack().fillna(0)

    # Rename the index and columns as needed
    # df_grouped.index = 'Cluster_' + df_grouped.index.astype(str)
    # df_grouped.columns = df_grouped.columns.map(lambda x: str(x).zfill(3))

    # Generate the boolean mask where all values in row are > 0
    valid_clusters = df_grouped.index[df_grouped.all(axis=1)]
    filtered_rows = []
    for cluster in valid_clusters:
        cluster_df = df_clstr[df_clstr['cluster'] == cluster]
        for org in cluster_df['Org'].unique():
            df_grouped = cluster_df[cluster_df['Org'] == org]
            max_identity_row = df_grouped.loc[df_grouped['identity'].idxmax()]
            filtered_rows.append(max_identity_row)
    
    filtered_df_clstr = pd.DataFrame(filtered_rows)
    filtered_df_clstr.reset_index(drop=True, inplace=True)
    filtered_df_clstr['cluster'] = 'Cluster_' + filtered_df_clstr['cluster'].astype(str).str.zfill(3)

             


    # Loop over each unique cluster
    for cluster in filtered_df_clstr['cluster'].unique():
        # Get the identifiers for this cluster
        identifiers = filtered_df_clstr[filtered_df_clstr['cluster'] == cluster]['identifier'].tolist()
        
        # Filter the fasta file based on these identifiers
        filter_fasta_file(paths["fasta_output"], identifiers, cluster, paths["core_clusters"])

    setup_future_paths(paths, rep_seq_path)

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="CD-HIT clustering")
#     parser.add_argument("-p","--path", help="Path where the cdhit program is located", default="/usr/bin")
#     parser.add_argument("-c", "--identity", type=float, default=0.9, help="Sequence identity threshold (default: 0.9)")
#     parser.add_argument("-d", "--length_diff", type=int, default=0, help="Length difference cutoff (default: 0)")
#     parser.add_argument("-sc", "--length_cutoff", type=int, default=1, help="Length difference cutoff in amino acid (default: 1)")
#     parser.add_argument("-aL", "--align_cov_long", type=float, default=0.9, help="Alignment coverage for longer sequence (default: 0.9)")
#     parser.add_argument("-aS", "--align_cov_short", type=float, default=0.9, help="Alignment coverage for shorter sequence (default: 0.9)")
#     parser.add_argument("-g", "--memory", type=int, default=1, help="Maximum available memory in GB (default: 1)")
#     parser.add_argument("-n", "--word_length", type=int, default=9, help="Word length (default: 3)")
#     args = parser.parse_args()
#     main(args)