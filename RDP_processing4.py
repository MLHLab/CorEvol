"""
Created on Fri May  5 14:39:01 2023

Parses RDP output in multiprocess functioning.

@author: abhishake
"""
import os
import glob
import pandas as pd
from openrdp import Scanner
import argparse
import json
from concurrent.futures import ProcessPoolExecutor

def path_setup(config_file = 'file_paths.json'):
    with open('file_paths.json', 'r') as f:
        paths = json.load(f)
    return paths

def setup_future_paths(paths, nonRdp, Rdp):
    paths["nonRdp_list"] = nonRdp
    paths["Rdp_list"] = Rdp
    with open('file_paths.json', 'w') as f:
        json.dump(paths, f)

def process_fasta_file(fasta_file, scanner_config, test_count):
    try:
        scanner = Scanner(cfg=scanner_config, methods=('bootscan', 'maxchi', 'siscan', 'chimaera', 'threeseq', 'rdp'))

        # Run scans on each .fa file
        results = scanner.run_scans(fasta_file)
        
        # Write the results to a temporary CSV file
        temp_csv = fasta_file.replace('.fa', '_results.csv')
        with open(temp_csv, 'w') as outfile:
            results.write(outfile)
        # print(fasta_file)

        # Read the results into a DataFrame
        df = pd.read_csv(temp_csv)
        # Filter the DataFrame based on Pvalue
        filtered_df = df[df['Pvalue'] < 0.05]
        
        # Group by 'Recombinant' and count the number of unique 'Method'
        grouped = filtered_df.groupby('Recombinant')['Method'].nunique()

        # Select the 'Recombinant' entries that appear for at least 4 different 'Method'
        result = grouped[grouped >= test_count].index
        
        # Check if the filtered DataFrame is empty
        if result.empty:
            # Return the filename without the path and extension as non-RDP
            return os.path.basename(fasta_file).replace('.fa', ''), False
        else:
            # Return the filename without the path and extension as RDP
            return os.path.basename(fasta_file).replace('.fa', ''), True
    except :
        # print(f"Error processing file: {fasta_file}")
        # print(f"Error message: {str(e)}")
        # print("Skipping to the next file...")
        return None, None

def process_fasta_file_with_retries(fasta_file, scanner_config, test_count, max_retries=3):
    for attempt in range(max_retries):
        try:
            return process_fasta_file(fasta_file, scanner_config, test_count)
        except Exception as e:
            print(f"Error processing file: {fasta_file}")
            print(f"Error message: {str(e)}")
            print(f"Attempt {attempt + 1} failed. Retrying...")
            if attempt == max_retries - 1:
                print("Max retries reached. Skipping to the next file...")
    return None, None

def main(rdp_config="/home/jisiasr/Abhishake/OpenRDP/tests/test_cfg_mod.ini", counts = 4, threads = 4):
    paths = path_setup('file_paths.json')

    scanner_config = rdp_config
    
    test_count = counts
    
    # Get a list of all .fa files in the /temp folder
    fasta_files = glob.glob(paths["aligned_core_cluster"] + "/*.fa")

    # Using ProcessPoolExecutor for parallel processing
    empty_results_files = []
    recom_result_files = []
    
    with ProcessPoolExecutor(max_workers=threads) as executor:
        results = list(executor.map(process_fasta_file_with_retries, fasta_files, [scanner_config] * len(fasta_files), [test_count] * len(fasta_files)))

    # Filter results to get files with no significant results
    for result, is_rdp in results:
        if result is not None:
            if is_rdp:
                recom_result_files.append(result)
            else:
                empty_results_files.append(result)

    # Write the list of filenames with empty result output to a file
    nonRdp = paths["temp_dir"] + "/nonRDP_list.txt"
    with open(nonRdp, 'w') as outfile:
        for filename in empty_results_files:
            outfile.write(f"{filename}\n")
    
    # Write the list of filenames with RDP to a file
    Rdp = paths["temp_dir"] + "/RDP_list.txt"
    with open(Rdp, 'w') as outfile:
        for filename in recom_result_files:
            outfile.write(f"{filename}\n")
    
    setup_future_paths(paths, nonRdp, Rdp)        



# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='Parallel processing of FASTA files.')
#     parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads (default: 4)')
#     parser.add_argument('-r','--rdp_config', default="/home/jisiasr/Abhishake/OpenRDP/tests/test_cfg_mod.ini", help='Path where internal parameters of RDP scanner is saved')
#     parser.add_argument('-x', '--counts', default=4, type=int, help='Number of different RDP testing methodology used to confidently conclude a sequence to be recombinant (default : 4, max : 6)')
#     args = parser.parse_args()
#     main(args)
    
    