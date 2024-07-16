#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 13:51:54 2024

@author: abhishake

The final pipeline function
"""

import argparse
import seq_input1
import cdhit_process2
import aligning3
import RDP_processing4
import CodeML_process52
import out_gen6

import subprocess

def run_script(script_path, *args):
    command = ['python', script_path] + list(args)
    subprocess.run(command, check=True)

def get_args():
    parser = argparse.ArgumentParser(description="Pipeline for running CorEvol.")
    parser.add_argument("-i", "--directory", required=True, help="Directory containing the FASTA files")
    parser.add_argument("-o", "--output", required=True, help="Directory containing the output files")
    parser.add_argument("-p", "--path", help="Path where the cdhit program is located", default="/usr/bin")
    parser.add_argument("-c", "--identity", type=float, default=0.9, help="Sequence identity threshold (default: 0.9)")
    parser.add_argument("-d", "--length_diff", type=int, default=0, help="Length difference cutoff (default: 0)")
    parser.add_argument("-sc", "--length_cutoff", type=int, default=1, help="Length difference cutoff in amino acid (default: 1)")
    parser.add_argument("-aL", "--align_cov_long", type=float, default=0.9, help="Alignment coverage for longer sequence (default: 0.9)")
    parser.add_argument("-aS", "--align_cov_short", type=float, default=0.9, help="Alignment coverage for shorter sequence (default: 0.9)")
    parser.add_argument("-g", "--memory", type=int, default=1, help="Maximum available memory in GB (default: 1)")
    parser.add_argument("-n", "--word_length", type=int, default=9, help="Word length (default: 3)")
    parser.add_argument("-pC", "--phylogeny_cutoff", type=int, default=300, help="Minimum length of nucleotides prior to alignment (default : 500)")
    parser.add_argument('-r', '--rdp_config', default="/home/jisiasr/Abhishake/OpenRDP/tests/test_cfg_mod.ini", help='Path where internal parameters of RDP scanner is saved')
    parser.add_argument('-x', '--counts', default=4, type=int, help='Number of different RDP testing methodology used to confidently conclude a sequence to be recombinant (default: 4, max: 6)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads (default: 4)')
    return parser.parse_args()


def pipeline(args):
    # Run dummy1Up
    seq_input1.main(directory=args.directory, result=args.output)

    # # Run dummy2Up
    cdhit_process2.main(path=args.path, identity=args.identity, length_diff=args.length_diff,
                  length_cutoff=args.length_cutoff, align_cov_long=args.align_cov_long,
                  align_cov_short=args.align_cov_short, memory=args.memory, word_length=args.word_length)

    # # Run dummy2_3
    aligning3.main(phylogeny_cutoff=args.phylogeny_cutoff)
    
    # # # # Run dummy3_Parallel
    RDP_processing4.main(threads=args.threads, rdp_config=args.rdp_config, counts=args.counts)
    
    # # # Run dummy4Up
    CodeML_process52.main()
    # # run_script('CodeML_process5.py')
    
    # # # Run dummy5
    out_gen6.main()
    

if __name__ == "__main__":
    args = get_args()
    pipeline(args)
