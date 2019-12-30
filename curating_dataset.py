# This script reads all the .txt files and converts into csv by reading first 3 columns.
# There are two new columns in the csv file for the metric and the target values.
#
# Author: Md Mahmudulla Hassan
# Date: 09/02/2018

import os
import glob
import pandas as pd

# Directory paths
output_dir = "dataset/actives/curated"
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

DATA_DIR = "dataset/actives/txt"
DATA_FILES = glob.glob(os.path.join(DATA_DIR, "*.txt"))

print("Total number of proteins: ", len(DATA_FILES))

# Curating the files
for _file in DATA_FILES:
    try:
      df = pd.read_csv(_file, sep='\t', header=None, usecols=[0, 1, 2], names = ["SMILES", "ChEMBL", "target"])
      # Drop rows if the target is None
      df = df.dropna(axis=0, subset=["target"])
      split_column = lambda x : pd.Series([i for i in x.split("=")])
      df['metric'] = df.target.apply(split_column)[0]
      df['target'] = df.target.apply(split_column)[1]
      dir_name, file_name = os.path.split(_file)
      df.to_csv(os.path.join(dir_name, output_dir, os.path.splitext(file_name)[0]+'.csv'), index=False)
    except Exception as e:
      print("ERROR IN " + _file)
      print(e)
      continue
