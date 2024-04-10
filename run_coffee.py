'''
Author: Musaddiq Lodi 
Date Created: 1/17/2023
Organization: Biological Networks Laboratory, VCU
Purpose: Partition scRNA-seq matrix into separate GRN inference; part of grn_clustering package
'''
import os
import pandas as pd 
import os 
import numpy as np
from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.cluster import SpectralClustering
import scanpy as sc 
import sys
import argparse


def runGRN(input_dir, home_dir, exp_file, coffee_threshold):

    #print(os.getcwd())

    print("Calculating Consensus GRN for each Cluster. This step may take a while...")
    os.chdir(home_dir)
    print("Calculating GRNs for", exp_file, "...")
    print("Creating Top 500 Gene File for file", exp_file)
    os.system(f"Rscript {input_dir}/grn_algs/tradeseq_ordering.r {os.getcwd()} {exp_file} > /dev/null 2>&1")
    print("Running PPCOR...")
    os.system(f"Rscript {input_dir}/grn_algs/ppcor_run.r {os.getcwd()} > /dev/null 2>&1")
    print("PPCOR Completed.")
    print("Running GENIE3...")
    os.system(f"Rscript {input_dir}/grn_algs/genie3_run.r {os.getcwd()} > /dev/null 2>&1")
    print("GENIE3 Completed.")
    print("Running GRNBOOST2...")
    os.system(f"python {input_dir}/grn_algs/grnboost2_run.py {os.getcwd()} > /dev/null 2>&1")
    print("Completed individual inference run for each GRN.")
    print("Running COFFEE consensus algorithm")
    #print(f"Python {input_dir}/coffee_consensus.py {os.getcwd()}/grn_networks 0.65 {os.getcwd()}/grn_networks")
    os.system(f"Python {input_dir}/grn_algs/coffee_consensus.py {os.getcwd()}/grn_networks {coffee_threshold} {os.getcwd()}/grn_networks > /dev/null 2>&1")
    print("COFFEE consensus completed")
    print(f"GRN inference for exp_file {exp_file} completed.")
        
        
def main():
    
    parser = argparse.ArgumentParser(description='Partition raw counts scRNA-seq expression into clusters and infer consensus GRN per cluster')
    
    # Add arguments
    parser.add_argument('--input_dir', default = os.getcwd(), help='Directory where scripts will be run from; default is the current working directory')
    parser.add_argument('--coffee_threshold', default = 0.65, help='Desired threshold for coffee algorithm, default is 0.65')
    parser.add_argument('--exp_file', help='Raw Counts scRNAseq file')
    parser.add_argument('--home_dir', help='Directory where output should be written to; default will be written to results within the input_dir')
    
    # Parse the command-line arguments
    args = parser.parse_args()

    # Access the arguments
    input_dir = args.input_dir
    coffee_threshold = args.coffee_threshold
    exp_file_path = args.exp_file
    home_dir = args.home_dir
    
    if not os.path.exists(home_dir):
        # If not, create the directory
        os.makedirs(home_dir)
        print(f"Directory '{home_dir}' created.")
    else:
        print(f"Directory '{home_dir}' already exists.")
    # Run the GRN inference
    runGRN(input_dir, home_dir, exp_file_path, coffee_threshold)
    
    
    print('\n')
    print("All GRN Inferences are complete.")
if __name__ == '__main__':
    main()