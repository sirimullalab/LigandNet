# Script for predicting Enamine 2M compounds using pre-extracted TPATF features
# This script does the following things:
#  1. Reads large csv files in chunks
#  2. For each of the chunks, create pre-defined number of processes and 
#  3. In each of the processes, reads the features and evaluates using all the models
#
# Author: Md Mahmudulla Hassan
# Department of Computer Science and School of Pharmacy, UTEP
# Last modifed: 10/16/2018

from __future__ import print_function
import pandas as pd
import numpy as np
import multiprocessing as mp
import os
from sklearn.externals import joblib
import sys

py_version = "27" if sys.version_info[0] < 3 else "35"
cores = 28 #mp.cpu_count()
partitions = cores
dataset_dir = "/data2/datasets/enamine"
fingerprint_files = ["Enamine_advance_TPATF_fingerprints.csv", "Enamine_hts_tpatf_fingerprint.csv"]
models_dir = "/data2/mhassan/LigandNet/Models/models_generated/all_models"
output_dir = "output/py" + py_version

def get_models():
    # Read the model names
    with open("/data2/mhassan/LigandNet/Models/models_generated/niners_py" + py_version + ".txt", 'r') as f:
        model_names = [line.split('\n')[0] for line in f.readlines()]
    names = [i.split("classifier_")[1] for i in model_names]
    # Yield models
    for name, model in zip(names, model_names):
        yield name, joblib.load(os.path.join(models_dir, model))


def work(data):
    # Predict SMILES in a chunk
    output = []
    data = data.apply(lambda x: x.tolist(), axis=1)
    ids = []
    smiles = []
    features = []
    for d in data:
        ids.append(d[0])
        smiles.append(d[1])
        features.append(np.array([float(i) for i in d[2].split(' ')], dtype=np.float32))

    features = np.array(features, dtype=np.float32)
    for model_name, model in get_models():
            # Predict and add to the output if there is a hit
            pred = []
            try:
                pred = model.predict_proba(features)
                pred = [i[1] for i in pred]
                hits = [i>0.5 for i in pred]
            except:
                pred = model.predict(features)
                hits = [i==1.0 for i in pred]

            if len(hits) > 0: output.extend([[model_name, i, s, p] for i, s, p, h in zip(ids, smiles, pred, hits) if h])

    return output

def parallelize(_file, chunk_number, data, func):
    data_split = np.array_split(data, partitions)
    pool = mp.Pool(cores)
    output = pool.map(func, data_split)
    pool.close()

    # with mp.Pool(cores) as pool:
    #     output = pool.map(func, data_split)
    
    # Filtering out empty results
    output = [i for i in output if len(i)>0]
    # Flatten the output
    result = []
    for i in output:
        result.extend(i)
    labels = ["protein", "id", "smiles", "prediction"]
    output_df = pd.DataFrame.from_records(result, columns=labels)
    output_file = os.path.join(output_dir, _file[:-4] + "_" + str(chunk_number) + ".csv")
    output_df.to_csv(output_file, index=None)

if __name__=="__main__":    
    for _file in fingerprint_files:
        data_file = os.path.join(dataset_dir, _file)
        print("Predicting ", data_file)
        for i, chunk in enumerate(pd.read_csv(data_file, chunksize=10**5)):
            parallelize(_file, i, chunk, work)
    
