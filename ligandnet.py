# Script for predicting ligand activity using LigandNet models
# Author: Md Mahmudulla Hassan
# Department of Computer Science and School of Pharmacy, UTEP
# Last modified: 02/20/2020

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
import argparse


class LigandNet(object):
    MODELS_DIR = os.path.join('models/files')

    def __init__(self):
        self.load_models()

    def load_models(self):
        # TODO: Avoid loading all the models
        # Read the best models
        with open('best_models.txt', 'r') as f:
            best_models = f.read().splitlines()

        self.uniprot_ids = [model_path[:6] for model_path in best_models]
        self.models = [joblib.load(os.path.join(
            self.MODELS_DIR, model_path)) for model_path in best_models]

    def get_features(self, input, input_type):
        # TODO: Add functionality for reading from a smi file containing a bulk of smiles
        ft = FeatureGenerator()
        if input_type == 'smiles':
            ft.load_smiles(input)
        else:
            ft.load_sdf(input)
        cmpd_id, features = ft.extract_tpatf()
        return cmpd_id, features.reshape(-1, 2692)

    # Get predictions
    def get_prediction(self, input, input_type, confidence_threshold=0.5):
        results = {}
        cmpd_id, features = self.get_features(input, input_type)
        cmpd_id = np.array(cmpd_id)
        for uniprot_id, model in tqdm(zip(self.uniprot_ids, self.models), total=703):
            # np.float32 is not json serializable, take float64
            pred = model.predict_proba(features).astype(np.float64)[:, 1].round(2)
            mask = pred >= confidence_threshold
            for _id, _pred in zip(cmpd_id[mask], pred[mask]):
                # Create a dictionary for compound if not exists
                if _id not in results.keys():
                    results[_id] = {}
                # Update the compound result dictionary
                results[_id].update({uniprot_id: _pred})
        return results


if __name__ == "__main__":
    start = time.time()
    parser = argparse.ArgumentParser(
        description="Ligand activity prediction using LigandNet")
    parser.add_argument('--sdf', action='store',
                        dest='sdf', help='SDF file location')
    parser.add_argument('--smiles', action='store', type=str,
                        dest='smiles', help='SMILES')
#     parser.add_argument('--out', action='store', dest='out',
#                         required=False, help='Output directory')
    parser.add_argument('--confidence', action='store', dest='confidence', type=float,
                        default=0.50, help='Minimum confidence to consider for prediction. Default is 0.5')

    args = parser.parse_args()

    if not (args.smiles or args.sdf):
        parser.error('No input found. Provide --smiles or --sdf')

    print(f"Loading the LigandNet models ...")
    l = LigandNet()

    if args.sdf is not None:
        if not os.path.isfile(args.sdf):
            raise FileNotFoundError(errno.ENOENT, os.strerror(
                errno.ENOENT), args.sdf)

        results = l.get_prediction(args.sdf, 'sdf', args.confidence)
        print(results)

    if args.smiles is not None:
        results = l.get_prediction(args.smiles, 'smiles', args.confidence)
        print(results)
