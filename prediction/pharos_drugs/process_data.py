import pandas as pd
import json
import sys
sys.path.append("../../../ddt/")
from utility import FeatureGenerator
import numpy as np
from tqdm import *
import multiprocessing as mp
from functools import partial

def get_features(row):
    smiles = row['smiles']
    ft = FeatureGenerator()
    ft.load_smiles(smiles)
    try:
        _, features = ft.extract_tpatf()
        return features
    except: return None


def parallelize(data, func, num_of_processes=mp.cpu_count()):
    data_split = np.array_split(data, num_of_processes)
    pool = mp.Pool(num_of_processes)
    data = pd.concat(pool.map(func, data_split))
    pool.close()
    pool.join()
    return data

def run_on_subset(func, data_subset):
    return data_subset.apply(func, axis=1)

def parallelize_on_rows(data, func, num_of_processes=mp.cpu_count()):
    return parallelize(data, partial(run_on_subset, func), num_of_processes)


pharos_drugs = pd.read_csv("pharos_drug_activity.csv", usecols=[0, 1, 6])
pharos_drugs.dropna(inplace=True)
pharos_drugs['features'] = parallelize_on_rows(pharos_drugs, get_features)
pharos_drugs.dropna(inplace=True)
f = np.array(pharos_drugs.features.apply(lambda x: x.flatten().tolist()).values.tolist())
np.save('features', f)
pharos_drugs['id', 'target_id', 'smiles'].to_csv("pharos_drugs.csv")
