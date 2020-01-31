from __future__ import print_function
import pandas as pd
import numpy as np
import os
from sklearn.externals import joblib
import sys
sys.path.append("../../../ddt/")
from utility import FeatureGenerator
from glob import glob
from sklearn.externals import joblib
import pickle
import pandas as pd
import errno
import json
from collections import OrderedDict
import time

MODELS_DIR = os.path.abspath("../../models")
MODELS_INFO = pd.read_csv("../../all_models_info.csv")
PY_VERSION = "py27" if sys.version_info[0] < 3 else "py35"

# Load the models
def get_models():
    for model_name, model_type in MODELS_INFO[MODELS_INFO['version'] == PY_VERSION][['name', 'model_type']].values:
        model_file = os.path.join(MODELS_DIR, model_name)
        if model_type == "RandomForestClassifier": model = joblib.load(model_file + '.rfc')
        elif model_type == "SVC": model = joblib.load(model_file + '.svc')
        elif model_type == "MLPClassifier": model = joblib.load(model_file + '.mlp')  
        elif model_type == 'XGBClassifier': model = pickle.load(open(model_file + '.xgb', 'rb'))
        yield model_name, model_type, model


# Read the zinc ids from sdf files
def get_zinc_id(_file):
    zinc_ids = []
    content = None
    with open(_file, 'r') as f:
        content = f.readlines()
    if content is None: return None
    content = [line[:-1] for line in content]
    return [line for line in content if line[:4] == 'ZINC']


if __name__=="__main__":
    start = time.time()
    import argparse
    parser = argparse.ArgumentParser(description="Prediction script for zinc database")    
    parser.add_argument('--sdf', action='store', dest='sdf', required=True, help='SDF file location')
    parser.add_argument('--out', action='store', dest='out', required=True, help='Output directory')
    parse_dict = vars(parser.parse_args())

    if not os.path.isfile(parse_dict['sdf']):
    	raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), parse_dict['sdf'])
    
    ft = FeatureGenerator()
    ft.load_sdf(parse_dict['sdf'])
    features = ft.extract_tpatf()
    features = np.array(features).reshape((-1, 2692))

    # predict using all the models
    results = OrderedDict()
    for _name, _type, _model in get_models():
        print("Predicting using {}".format(_name))
        _, sdf_filename = os.path.split(parse_dict['sdf'])
        
        results['zinc_ids'] = get_zinc_id(parse_dict['sdf'])
        results['compound_count'] = features.shape[0]

        predict_info = OrderedDict()
        predict_info['model_type'] = _type
        _pred = _model.predict(features) if _type=='SVC' else _model.predict_proba(features)[:, 1]
        predict_info['activity'] = _pred.tolist()
        results[_name] = predict_info

        # Write the output
        with open(os.path.join(parse_dict['out'], sdf_filename + ".csv"), 'w') as f:
            json.dump(results, f)
        #break
    end = time.time()
    print("Execution time: {} seconds.".format(round(end - start, 2)))