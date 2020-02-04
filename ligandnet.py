import os
import joblib
import argparse
import time
import pandas as pd
import numpy as np
from ddt.utility import FeatureGenerator
import errno
import json
from collections import OrderedDict
from tqdm import tqdm

MODELS_DIR = os.path.join('models/files')


# Load the models
def get_models():
    # Read the best models
    with open('best_models.txt', 'r') as f:
        best_models = f.read().splitlines()

    for model_path in best_models:
        # if model_path[-3:] == 'svc': continue
        yield [model_path[:6], joblib.load(os.path.join(MODELS_DIR, model_path))]


def get_features(_input, input_type):
    ft = FeatureGenerator()
    if input_type == 'smiles':
        ft.load_smiles(_input)
    else:
        ft.load_sdf(_input)
    _id, features = ft.extract_tpatf()
    return zip(_id, features)


# Get predictions
def get_prediction(_input, input_type, confidence_threshold=0.5):
    # predict using all the models
    results = {}
    for (_id, features) in get_features(_input, input_type):
        predict_info = {}
        for [_name, _model] in tqdm(get_models(), total=703):
            # print("Finding ligand activity for {}".format(_name))
            _pred = _model.predict_proba(features.reshape(1, -1))[:, 1]
            if _pred[0] < confidence_threshold:
                continue
            predict_info[_name] = _pred[0]

        if len(predict_info) > 0:
            results[_id] = predict_info 

    print(results)
    return results


if __name__ == "__main__":
    import argparse
    start = time.time()
    import argparse
    parser = argparse.ArgumentParser(
        description="Ligand activity prediction using LigandNet")
    parser.add_argument('--sdf', action='store',
                        dest='sdf', help='SDF file location')
    parser.add_argument('--smiles', action='store',
                        dest='smiles', help='SMILES')
    parser.add_argument('--out', action='store', dest='out',
                        required=False, help='Output directory')
    parser.add_argument('--confidence', action='store', dest='confidence', type=float,
                        default=0.50, help='Minimum confidence to consider for prediction')
    parse_dict = vars(parser.parse_args())

    if parse_dict['sdf'] is not None:
        if not os.path.isfile(parse_dict['sdf']):
            raise FileNotFoundError(errno.ENOENT, os.strerror(
                errno.ENOENT), parse_dict['sdf'])

        get_prediction(parse_dict['sdf'], 'sdf', parse_dict['confidence'])

    if parse_dict['smiles'] is not None:
        get_prediction(parse_dict['smiles'],'smiles', parse_dict['confidence'])
