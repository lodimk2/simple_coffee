# Run GRNBOOST2 GRN Inference Algorithm

import sys
import pandas as pd
import os 
import scanpy as sc
from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names
from distributed import Client, LocalCluster
client = Client(processes = False)    

input_dir = sys.argv[1]

os.chdir(input_dir)

inDF = pd.read_csv("cluster_top_exp.csv", index_col= 0)

print(inDF.head())

adata = sc.AnnData(X=inDF.T)

sc.pp.log1p(adata)

# Get the log-transformed matrix as a Pandas DataFrame
log_transformed_df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
print("Printing log-transformed-df")
print(log_transformed_df.head())
print("Printing log-transformed-df to numpy")
print(log_transformed_df.to_numpy())

print(log_transformed_df.to_numpy().shape)

network = grnboost2(log_transformed_df.to_numpy(), client_or_address = client, gene_names = inDF.index)
network.columns = ['Gene1', 'Gene2', 'EdgeWeight']
print(network.head())


# Check if the directory exists
if not os.path.exists("grn_networks"):
    # If not, create the directory
    os.makedirs("grn_networks")
    print(f"Directory grn_networks created.")
else:
    print(f"Directory grn_networks already exists.")

network.to_csv("grn_networks/grnboost_network.csv", sep = "\t", index = False)
