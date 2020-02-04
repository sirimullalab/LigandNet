import json
import sys
import numpy as np
from tqdm import *
import random
import pandas as pd
import numpy as np
import sys
from sklearn.externals import joblib
import os
import json
import warnings
import time
warnings.filterwarnings("ignore")

home_dir = '../..'
MODELS_DIR = '../../uniprot_models'
PY_VERSION = '27' if sys.version_info[0] < 3 else '35'


# Load the models
def get_models():    
    model_names = open(os.path.join(home_dir, 'py' + PY_VERSION + '_uniprot_models.txt'), 'r').readlines()[0].split(',')
    for model in model_names:
        with open(os.path.join(MODELS_DIR, model), 'rb') as f:
            yield model, joblib.load(f) #, mmap_mode='r+')


def get_prediction(features):
    confidence = 0.9
    actives = []
    for model_name, model in get_models():
        if type(model).__name__ == "SVC": 
            pred = model.predict(features)
            pred = pred.reshape((-1, 1))
            if pred[:, 0] > confidence: actives.append(model_name)
        else:
            pred = model.predict_proba(features)
            pred = pred.reshape((-1, 2))
            if pred[:, 1] > confidence: actives.append(model_name)
            
    return actives

print("READING DATA...")
pharos_drugs = pd.read_csv("pharos_drugs.csv")
print("LOADING FEATURES...")
features = np.load('features.npy')
print("RUNNING MODELS...")
predictions = []
for feat in tqdm(features):
    #TODO: don't add list, add string
    predictions.append(get_prediction(feat[np.newaxis, :]))
print("SAVING RESULTS...")
pharos_drugs['actives'] = pd.Series.from_array(predictions)
pharos_drugs.to_csv('pharos_drugs_actives_'+PY_VERSION+'.csv', index=None)
