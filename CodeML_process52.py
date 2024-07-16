#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 18:49:16 2024

@author: abhishake

Code for codeml calculation
"""
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import numpy as np
import pandas as pd
# import argparse
import json
from Bio.Phylo.PAML import codeml
import multiprocessing as mp
import os
import signal

def path_setup(config_file = 'file_paths.json'):
    with open('file_paths.json', 'r') as f:
        paths = json.load(f)
    return paths
def setup_future_paths(paths, summary_file, tree_out):
    paths["summary_file"] = summary_file
    paths["concate_tree"] = tree_out
    with open('file_paths.json', 'w') as f:
        json.dump(paths, f)

def set_up_codeml(alignment_file, tree_file, working_dir):
    """
    Sets up the Codeml instance with specified options.
    
    Parameters:
        alignment_file (str): Path to the alignment file.
        tree_file (str): Path to the tree file.
        working_dir (str): Working directory for Codeml.
    
    Returns:
        Codeml: Configured Codeml instance.
    """
    # from Bio.Phylo.PAML import codeml

    cml = codeml.Codeml()
    cml.alignment = alignment_file
    cml.tree = tree_file
    cml.out_file = "results.txt"
    cml.working_dir = working_dir
    cml.set_options(seqtype=1, verbose=1, runmode=-2, cleandata=1, model=0,
                    NSsites=[0], fix_kappa=0, kappa=2, fix_omega=0, omega=1,
                    clock=1, getSE=1, RateAncestor=0, method=0, fix_alpha=1,
                    alpha=0.0, Malpha=0, ncatG=4, CodonFreq=2)
    return cml

def extract_dN_dS(results):
    """
    Extracts dN and dS values from the results of a Codeml run.
    
    Parameters:
        results (dict): The results dictionary returned from Codeml run.
    
    Returns:
        tuple: Two lists containing dN and dS values respectively.
    """
    dN = []
    dS = []
    if 'pairwise' in results:
        pairwise = results.get("pairwise")
        for key1, sub_dir1 in pairwise.items():
            for key2, sub_dir2 in sub_dir1.items():
                dn_value = sub_dir2.get("dN")
                ds_value = sub_dir2.get("dS")
                if dn_value is not None:
                    dN.append(dn_value)
                if ds_value is not None:
                    dS.append(ds_value)
    return dN, dS


def calculate_stats(data):
    """
    Calculates the mean and standard error of the given data array.

    Parameters:
        data (np.ndarray): An array of data points (floats or integers).

    Returns:
        tuple: A tuple containing the mean and standard error of the data.
    """
    mean_value = np.mean(data)
    # Standard Error (SE) = Standard Deviation (SD) / sqrt(number of samples)
    standard_error = np.std(data, ddof=1) / np.sqrt(len(data))
    return round(mean_value,5), round(standard_error,5)




def process_fasta_file(fasta_file):
    """
    Processes a FASTA file to format it for phylogenetic analysis, constructs a phylogenetic tree,
    and calculates non-synonymous (dN) and synonymous (dS) substitution rates.

    This function performs several operations:
    1. Reads an alignment from a FASTA file.
    2. Writes the alignment to a new file in PHYLIP format.
    3. Calculates pairwise distances between sequences and uses these to construct a neighbor-joining tree.
    4. Writes the tree to a file in Newick format.
    5. Sets up and runs a CodeML analysis using the PHYLIP file and the tree file, extracting the dN and dS values.

    Parameters:
        fasta_file (str): The path to the FASTA file to be processed.

    Returns:
        tuple: A tuple containing two lists, dN and dS, which are the non-synonymous and synonymous
               substitution rates respectively.

    Raises:
        IOError: If the file cannot be opened, read, or written.
        ValueError: If the input data is improperly formatted or insufficient for analysis.

    Example:
        >>> dN, dS = process_fasta_file("path/to/fasta_file.fa")
        >>> print(f"dN: {dN}, dS: {dS}")
    """
    # Read the alignment from fasta
    try :
        alignment = AlignIO.read(fasta_file, "fasta")

        # Efficient file write
        output_file_name = fasta_file.replace(".fa", "_formatted.phy")
        with open(output_file_name, "w") as output_file:
            num_sequences = len(alignment)
            alignment_length = alignment.get_alignment_length()
            output_file.write(f"{num_sequences} {alignment_length}\n")
            for record in alignment:
                sequence_name = record.id[:10].ljust(10)
                sequence_data = f"{sequence_name}    {record.seq}\n"
                output_file.write(sequence_data)

        # Calculate distances and construct tree
        calculator = DistanceCalculator('identity')
        constructor = DistanceTreeConstructor(calculator, 'nj')
        tree = constructor.build_tree(alignment)
        tree.ladderize()

        # Efficient tree writing
        tree_file_name = fasta_file.replace(".fa", "_tree.nwk")
        with open(tree_file_name, "w") as file:
            Phylo.write(tree, file, "newick")

        # Setting up CodeML
        cml = set_up_codeml(output_file_name, tree_file_name, ".")
        results = cml.run()
        dN, dS = extract_dN_dS(results)

        return dN, dS
    except Exception as e:
        print(f"Error processing file {fasta_file}: {str(e)}")
        return [], []  # Return empty lists for dN and dS
    

def write_results_to_tsv(output_filename, results):
    with open(output_filename, 'w') as file:
        file.write("Filename\tdN Average\tSE of dN\tdS Average\tSE of dS\tdN/dS\n")
        for result in results:
            file.write(f"{result['filename']}\t{result['dN_avg']}\t{result['dN_se']}\t{result['dS_avg']}\t{result['dS_se']}\t{result['ratio']}\n")


def retry_process_file(fasta_file, max_retries=3):
    """
    Attempts to process the provided FASTA file up to a specified number of retries.
    
    Parameters:
        fasta_file (str): Path to the FASTA file.
        max_retries (int): Maximum number of retries.
    
    Returns:
        tuple: dN and dS lists if successful, otherwise empty lists.
    """
    attempts = 0
    while attempts < max_retries:
        dN, dS = process_fasta_file(fasta_file)
        # Check if results are valid (not empty and not containing only NaNs)
        if dN and dS and not (np.isnan(dN).all() or np.isnan(dS).all()):
            return dN, dS
        attempts += 1
        print(f"Retry {attempts}/{max_retries} for file {fasta_file}.")
    return [], []  # Return empty lists if all retries fail

def process_fasta_file_wrapper(fasta_file):
    print(f"Processing file: {fasta_file}")
    
    # Set a timeout of 5 minutes (adjust as needed)
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(300)  # 5 minutes
    
    try:
        # dN, dS = retry_process_file(fasta_file)
        dN, dS = process_fasta_file(fasta_file)
        
        print("Processing completed. Calculating stats...")
        dN_avg, dN_se = calculate_stats(dN)
        dS_avg, dS_se = calculate_stats(dS)
        dN_dS_ratio = dN_avg / dS_avg if dS_avg != 0 else float('inf')
        
        print(f"Finished processing file: {fasta_file}")
        
        return {
            'filename': fasta_file,
            'dN_avg': dN_avg,
            'dN_se': dN_se,
            'dS_avg': dS_avg,
            'dS_se': dS_se,
            'ratio': dN_dS_ratio
        }
    except TimeoutError:
        print(f"Timeout occurred for file: {fasta_file}")
        return {
            'filename': fasta_file,
            'dN_avg': None,
            'dN_se': None,
            'dS_avg': None,
            'dS_se': None,
            'ratio': None
        }
    finally:
        signal.alarm(0)  # Cancel the alarm
        
def timeout_handler(signum, frame):
    raise TimeoutError("Function execution timed out.")

def concatenate_aligned_fasta(file_list, output_file, header_map):
    # Dictionary to store the concatenated sequences for each organism
    concatenated_seqs = {}

    # Iterate over each file in the list
    for file_name in file_list:
        # Parse the aligned FASTA file using Biopython's SeqIO
        for record in SeqIO.parse(file_name, 'fasta'):
            # Split the sequence header to get the organism name
            org_name = record.id.split("_")[0]
            seq_data = str(record.seq)

            # Add the sequence data to the dictionary for the corresponding organism
            concatenated_seqs.setdefault(org_name, '')
            concatenated_seqs[org_name] += seq_data

    # Create a list to store the concatenated sequence records
    concatenated_records = []

    # Create SeqRecord objects for each concatenated sequence
    for org_name, seq_data in concatenated_seqs.items():
        # Check if the organism name exists in the header map
        if org_name in header_map:
            new_header = header_map[org_name]
        else:
            new_header = org_name

        concatenated_record = SeqRecord(Seq(seq_data), id=new_header, description='')
        concatenated_records.append(concatenated_record)

    # Write the concatenated sequences to the output file
    with open(output_file, 'w') as output:
        SeqIO.write(concatenated_records, output, 'fasta')

    print(f"Concatenated aligned FASTA file created: {output_file}")

   



def main(threads=4):
    paths = path_setup('file_paths.json')
    
    # Define the path where the FASTA files are located
    path_prefix = paths["aligned_core_cluster"]

    # Read the FASTA file names from the text file
    fasta_files = []
    with open(paths["nonRdp_list"], "r") as file:
        for line in file:
            # Strip newline and any surrounding whitespace
            trimmed_name = line.strip()
            # Construct the full file path and add it to the list
            if not trimmed_name.endswith(".fa"):
                trimmed_name += ".fa"
            full_path = path_prefix + "/" + trimmed_name
            fasta_files.append(full_path)
            
    filename_map = pd.read_csv(paths['filename_map'], sep = "\t", dtype=str, header=0)
    filename_map['original'] = filename_map['original'].apply(lambda x: os.path.splitext(x)[0])
    org_map = dict(zip(filename_map['coded'], filename_map['original']))
    tree_out = paths["result_dir"] + "/concatenate.tree"
    concatenate_aligned_fasta(fasta_files,tree_out,org_map)

    # Now fasta_files list contains the full paths with the adjusted names
    # print("FASTA files to process:", fasta_files)

    results = []

    for fasta_file in fasta_files:
        dN, dS = retry_process_file(fasta_file)
        dN_avg, dN_se = calculate_stats(dN)
        dS_avg, dS_se = calculate_stats(dS)
        dN_dS_ratio = dN_avg / dS_avg if dS_avg != 0 else float('inf')  # Avoid division by zero

        results.append({
            'filename': fasta_file,
            'dN_avg': dN_avg,
            'dN_se': dN_se,
            'dS_avg': dS_avg,
            'dS_se': dS_se,
            'ratio': dN_dS_ratio
        })
        
    # with mp.Pool(processes=1) as pool:
    #     results = pool.map(process_fasta_file_wrapper, fasta_files)
    
    summary_file = paths["temp_dir"] + "/result_summary_R1.tsv"

    write_results_to_tsv(summary_file, results)
    
    
    # print("Results written to 'results_summary.tsv'")
    setup_future_paths(paths, summary_file, tree_out)

# if __name__=='__main__':
#     parser = argparse.ArgumentParser(description='codeml implementation of core, non-redundant, non-recombinant sequences.')
#     args = parser.parse_args()
#     main(threads=4)
