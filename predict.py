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

MODELS_DIR = os.path.join('models/files')


# Get the model types
def get_model_name(ext):
    name_dict = {'rf': 'RandomForestClassifier',
                 'svc': 'SupportVectorClassifier', 'mlp': 'MLPClassifier', 'xgb': 'XGBoost'}
    return name_dict[ext]


# Load the models
def get_models():
    # Read the best models
    with open('best_models.txt', 'r') as f:
        best_models = f.read().splitlines()

    for model_path in best_models:
        # if model_path[-3:] == 'svc': continue
        yield [model_path[:6], joblib.load(os.path.join(MODELS_DIR, model_path))]


# Get predictions
def get_prediction(sdf_file, output=None):
    ft = FeatureGenerator()
    ft.load_sdf(sdf_file)
    zinc_ids, features = ft.extract_tpatf()

    # predict using all the models
    results = OrderedDict()
    for [_name, _model] in get_models():
        print("Predicting using {}".format(_name))
        results['zinc_ids'] = zinc_ids
        results['compound_count'] = len(zinc_ids)

        predict_info = OrderedDict()
        _pred = _model.predict_proba(features)[:, 1]
        predict_info['activity'] = _pred.tolist()
        results[_name] = predict_info

    if output is not None:
        if not os.path.isdir(output):
            os.makedirs(output)
        # Write the output
        _, sdf_file = os.path.split(sdf_file)
        sdf_file, _ = os.path.splitext(sdf_file)
        output_file = os.path.join(output, sdf_file + ".json")
        with open(output_file, 'w') as f:
            json.dump(results, f)
        print(f"Results are saved to {output_file}")

    return results


if __name__ == "__main__":
    import argparse
    start = time.time()
    import argparse
    parser = argparse.ArgumentParser(
        description="Prediction script for zinc database using LigandNet models")
    parser.add_argument('--sdf', action='store', dest='sdf',
                        required=True, help='SDF file location')
    parser.add_argument('--out', action='store', dest='out',
                        required=True, help='Output directory')
    parse_dict = vars(parser.parse_args())

    if not os.path.isfile(parse_dict['sdf']):
        raise FileNotFoundError(errno.ENOENT, os.strerror(
            errno.ENOENT), parse_dict['sdf'])

    # Predict
    results = get_prediction(parse_dict['sdf'], parse_dict['out'])
    end = time.time()
    print("Execution time: {} seconds.".format(round(end - start, 2)))
