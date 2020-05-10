import numpy as np
from sklearn.externals import joblib
import os
import json

MODELS_DIR = "../../models"
import sys

sys.path.append("../../../ddt/")
from utility import FeatureGenerator

PY_VERSION = "27" if sys.version_info[0] < 3 else "35"

# Extract features
ft = FeatureGenerator()
ft.load_sdf("20190509-L1300-FDA-approved-Drug-Library.sdf")
# ft.load_sdf("selleck_sample.sdf")
compound_list, features = ft.extract_tpatf()
if len(compound_list) != features.shape[0]:
    print("ERROR: Mismatch in total number of compounds and list of features")
    print(
        "       Number of compounds in the list: {} and feature cound: {}".format(
            len(compound_list), features.shape[0]
        )
    )
    sys.exit(1)

# features = np.load('drug_central.npy')
# compound_list = open('compound_list.txt').readlines()[0].split(',')


def get_models():
    model_names = (
        open("../../py" + PY_VERSION + "_models.txt", "r").readlines()[0].split(",")
    )
    for model in model_names:
        yield model, joblib.load(os.path.join(MODELS_DIR, model))


def get_prediction(confidence=0.9):
    results = {}
    for model_name, model in get_models():
        print("Predicting using {}".format(model_name))
        if type(model).__name__ == "SVC":
            pred = model.predict(features)
            mask = np.where(pred)
            if np.any(mask):
                results[model_name] = np.array(compound_list)[mask].tolist()
        else:
            pred = model.predict_proba(features)
            mask = np.where(pred[:, 1] > confidence)
            if np.any(mask):
                c_list = np.array(compound_list)[mask].tolist()
                results[model_name] = [
                    (c, str(f)) for c, f in zip(c_list, pred[:, 1][mask])
                ]
    return results


if __name__ == "__main__":
    results_dict = get_prediction()
    json.dump(results_dict, open("py" + PY_VERSION + "_results.json", "w"))
