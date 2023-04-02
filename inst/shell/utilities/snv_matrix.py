#!/usr/bin/python3

import os
import sys
from os import listdir
from os.path import isfile, join
import pandas as pd
import functools 
from collections import defaultdict

input_folder = sys.argv[1] # --input-dir
output_folder = sys.argv[2] # --output-dir

suffix = "_filtered_pass.vcf"
onlyfiles = [f for f in listdir(input_folder) if isfile(join(input_folder, f)) if f.endswith(suffix)]

min_shared_snv = 10
snv_dict = defaultdict(list)
sample_list = [] # to store sample names which will be column name of the dataframe
for f in onlyfiles:
    path_to_file = os.path.join(input_folder,f)
    current_file = open(path_to_file, 'r')
    sample_name = f.split("_")[0] 
    sample_list.append(sample_name)
    Lines = current_file.readlines()
    for line in Lines:
        if not line.startswith("#"):
            chrom = line.split()[0]
            position = line.split()[1]
            snv_index = chrom + "__" + position
            snv_dict[snv_index].append(sample_name)
            current_file.close()

snv_filtered = defaultdict(list)
 
for snv in snv_dict:
    if len(snv_dict[snv]) >= min_shared_snv:
        snv_filtered[snv] = snv_dict[snv]

snv_filtered_name_list = list(snv_filtered.keys())
snv_df = pd.DataFrame(0, index=snv_filtered_name_list, columns=sample_list)

for snv in snv_filtered:
    for sample in snv_filtered[snv]:
        snv_df.loc[snv, sample] = 1

snv_df.index.name = "SNVid"

path_to_saved_snv_df = join(output_folder, "SNV_matrix.csv")
snv_df.to_csv(path_to_saved_snv_df)
