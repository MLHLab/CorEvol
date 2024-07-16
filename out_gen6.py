#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 16:30:10 2024

@author: abhishake
"""
import pandas as pd
import numpy as np
import os
# import shutil
import json
from Bio import SeqIO

def path_setup(config_file = 'file_paths.json'):
    with open('file_paths.json', 'r') as f:
        paths = json.load(f)
    return paths


def actual_seq_name(df, rep_seq_list, header_df):
    # Create a dictionary to map coded sequence names to actual sequence names
    seq_name_map = dict(zip(header_df.iloc[:, 1], header_df.iloc[:, 0]))
    
    # Create a dictionary to map cluster IDs to coded sequence names
    cluster_seq_map = dict(zip(rep_seq_list['Cluster_ID'], rep_seq_list['identifier']))
    
    # Extract the cluster number from the 'Filename' column
    df['cluster_num'] = df['Filename'].str.extract(r'Cluster_(\d+)', expand=False)
    
    # Map cluster numbers to cluster IDs
    df['cluster_id'] = 'Cluster_' + df['cluster_num'].astype(str).str.zfill(3)
    
    # Map cluster IDs to coded sequence names
    df['coded_seq_name'] = df['cluster_id'].map(cluster_seq_map)
    
    # Map coded sequence names to actual sequence names
    df['Rep_seq_name'] = df['coded_seq_name'].map(seq_name_map)
    
    # Drop temporary columns
    df.drop(['cluster_num', 'cluster_id', 'coded_seq_name'], axis=1, inplace=True)
    
    return df

def replace_fasta_headers(fasta_file, header_df, output_dir):
    # Create a dictionary to map coded sequence names to actual sequence names
    seq_name_map = dict(zip(header_df.iloc[:, 1], header_df.iloc[:, 0]))
    
    # Parse the input FASTA file using SeqIO
    records = list(SeqIO.parse(fasta_file, 'fasta'))
    
    # Replace the coded headers with actual sequence names
    for record in records:
        coded_header = record.id
        actual_header = seq_name_map.get(coded_header, coded_header)
        record.id = actual_header
        record.description = ''
    output_file = output_dir + os.path.basename(fasta_file)
    # Write the modified records to the output FASTA file
    SeqIO.write(records, output_file, 'fasta')

def main():
    paths = path_setup('file_paths.json')
    
    # Step 1: Read the TSV file
    df = pd.read_csv(paths['summary_file'], sep='\t')
    rep_seq_list = pd.read_csv(paths['rep_seq_list'], sep = '\t')
    header_df = pd.read_csv(paths['csv_output'], names=["actual_seq_name","coded_seq_name"])
    # Step 2: Clean data by removing NaN or Inf values
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    
    # Modify the first column to extract just the cluster name part
    df['Filename'] = df['Filename'].apply(lambda x: x.split('/')[-1].split('.')[0][8:])
    
    df_mod = actual_seq_name(df,rep_seq_list,header_df)
    
    # Print the modified dataframe
    df_mod.to_csv("Selection_table.tsv", sep="\t", index=False)
    
    
    
    # setup_future_paths(paths, df_mod)
    
    # Step 3: Create directories for file classification
    positive_selection_dir = os.path.join(paths["result_dir"],'Positive_Selection/')
    purifying_selection_dir = os.path.join(paths["result_dir"],'Purifying_Selection/')
    neutral_selection_dir = os.path.join(paths["result_dir"],'Neutral_Selection/')
    recombination_dir = os.path.join(paths["result_dir"],'Recombination_Clusters/')
    os.makedirs(positive_selection_dir)
    os.makedirs(purifying_selection_dir)
    os.makedirs(neutral_selection_dir)
    os.makedirs(recombination_dir)
    
        
    # Step 4: Classify and copy files
    for index, row in df.iterrows():
        cluster_name = row['Filename']
        # Construct the original file path dynamically and check if it exists
        original_file_path = os.path.join(paths['core_clusters'], f'{cluster_name}.fa')
        
        if not os.path.exists(original_file_path):
            print(f"Warning: File not found {original_file_path}")
            continue
        
        dN_dS_ratio = row['dN/dS']
        
        if dN_dS_ratio > 1:
            # Files with positive selection
            replace_fasta_headers(original_file_path, header_df, positive_selection_dir)
            
            # shutil.copy(original_file_path, positive_selection_dir)
        elif (0 < dN_dS_ratio < 1):
            # Files with purifying selection
            replace_fasta_headers(original_file_path, header_df, purifying_selection_dir)
            # shutil.copy(original_file_path, purifying_selection_dir)
        elif (dN_dS_ratio == 0):
            # Files with purifying selection
            replace_fasta_headers(original_file_path, header_df, neutral_selection_dir)
    
    print("Files have been classified and copied to their respective folders.")
    
    #Step 5: copy RDP files
    # Define the path where the FASTA files are located
    path_prefix = paths["core_clusters"]

    # Read the FASTA file names from the text file
    # fasta_files = []
    with open(paths["Rdp_list"], "r") as file:
        for line in file:
            # Strip newline and any surrounding whitespace
            trimmed_name = line.strip()
            # Construct the full file path and add it to the list
            if not trimmed_name.endswith(".fa"):
                trimmed_name += ".fa"
            full_path = path_prefix + "/" + trimmed_name
            replace_fasta_headers(full_path, header_df, recombination_dir)
            

# if __name__ == "__main__":
#     main()
