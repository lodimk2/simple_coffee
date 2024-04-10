import pandas as pd
import os
import numpy as np
import sys
import shutil
from collections import defaultdict, Counter
directory = sys.argv[1]
threshold = sys.argv[2]
output_dir = sys.argv[3]

threshold_path = str(threshold) + "_consensus_net"
if os.path.exists(threshold_path) and os.path.isdir(threshold_path):
        shutil.rmtree(threshold_path)
def borda(ranked_lists):
    # Create a defaultdict to store the scores for each candidate
    edge_scores = defaultdict(int)
    
    # Count the number of candidates
    num_edges = len(ranked_lists[0])

    # Calculate the Borda scores using Counter
    for alg in ranked_lists:
        for i, edge in enumerate(alg):
            edge_scores[edge] += num_edges - i - 1

    # Sort the scores in descending order
    sorted_scores = dict(sorted(edge_scores.items(), key=lambda x: x[1], reverse=True))

    # Return the sorted scores
    return sorted_scores

def make_ranked_list(directory):

    file_list = []
    counter = 0
    min_thresh = 0
    for filename in os.listdir(directory):
        print("Going through filename", filename)
        if filename == ".DS_Store":
            continue
        if (".csv" not in filename) and (".tsv" not in filename):
            continue
        else:
            print("Got to this step!!")
            lct_file = os.path.join(directory, filename)
            print("lct_file", lct_file)
            lct_test_csv = pd.read_csv(lct_file, sep = "\t")
            lct_test_csv.columns = ['Gene1', 'Gene2', 'EdgeWeight']
            lct_test_csv = lct_test_csv[lct_test_csv['Gene1'] != lct_test_csv['Gene2']]
            print(lct_test_csv.head(5))
            line_count = len(lct_test_csv.index)
            if line_count > min_thresh:
                min_thresh = line_count
    for filename in os.listdir(directory):
        if filename == ".DS_Store":
            continue
        elif (".csv" not in filename) and (".tsv" not in filename):
            continue
        else:
            threshold = min_thresh
            #print("USING MIN THRESHOLD OF", threshold)
            f = os.path.join(directory, filename)
            test_csv = pd.read_csv(f, sep = "\t")
            test_csv.columns = ['Gene1', 'Gene2', 'EdgeWeight']
            test_csv = test_csv[test_csv['Gene1'] != test_csv['Gene2']]
            #print(test_csv)
            edges_1 = test_csv['Gene1'].tolist()[0:threshold]
            edges_2 = test_csv['Gene2'].tolist()[0:threshold]
            #weights = test_csv['EdgeWeight'].tolist()[0:threshold]
            edge_tuple = []
            edge_position_dict = {}

            for edge in range(len(edges_1)):
                edge_tuple.append((edges_1[edge], edges_2[edge]))

            file_list.append(edge_tuple)
            counter += 1
    return(file_list)


def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))



ranked_list = make_ranked_list(directory)


to_normalize = [i for i in range(len(ranked_list[0]))]
norm_list = list(NormalizeData(to_normalize))
norm_list.reverse()

result = borda(ranked_list)

edge_1 = []
edge_2 = []
scores = norm_list
for key in result.keys():
    edge_1.append(key[0])
    edge_2.append(key[1])

out_dict = {
    "Edge1" : edge_1,
    "Edge2" : edge_2, 
    "Scores" : scores
}

out_df = pd.DataFrame(out_dict)

# Write to csv
out_df = out_df[out_df.iloc[:, 2] >= float(threshold)]

os.chdir(output_dir)
threshold_path = str(threshold) + "_consensus_net"
if os.path.exists(threshold_path) and os.path.isdir(threshold_path):
        shutil.rmtree(threshold_path)
os.mkdir(threshold_path)
print("Created threshold_path of", threshold_path)
out_df.columns = ['Gene1', 'Gene2', 'EdgeWeight']
out_df.to_csv(threshold_path + "/consensus_network.tsv", sep = "\t", index = False)
print("Wrote final outfile")

