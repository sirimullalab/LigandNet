#################################################################
# This script trains Ligandnet models                           #
# Models are: XGBoost, RandomForest, Support Vector Classifier  #
#  and Neural Network                                           #
#                                                               #
# Author: Mahmudulla Hassan, UTEP CS and School of Pharmacy     #
# Last modified: 12/30/2019                                     #
#################################################################

#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import os
from glob import glob
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import metrics
import xgboost as xgb
import multiprocessing as mp
import pickle
import json
import os
import json
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import scorer, make_scorer
import joblib
import psutil
import numpy as np

classifier_loglevel = 0
gridsearch_loglevel = 2

# List the directories and files
active_dir = "pharos_database/actives_fingerprints"
decoy_dir = "pharos_database/decoys_fingerprints"

# Load uniprot ids
protein_to_uniprot = json.load(open('protein_to_uniprot.json', 'r'))
        
# active_files = glob(os.path.join(active_dir, '*.csv'))
# decoy_files = glob(os.path.join(decoy_dir, '*.csv'))
# print("Number of active and decoy files: {}, {}".format(len(active_files), len(decoy_files)))

class Train(object):
    
    def __init__(self, active_file, decoy_file, test_split=0.2, random_state=1, sample_threshold=20):
        self.active_file = active_file
        self.decoy_file = decoy_file
        self.test_split = test_split
        self.random_state = random_state
        self.sample_threshold = sample_threshold
        
        # Initialize directories
        self.model_dir = 'models'
        self.output_dir = 'reports'
        for _dir in [self.model_dir, self.output_dir]:
            if not os.path.isdir(_dir): os.makedirs(_dir)
        
        # Get the protein name and the uniprot id
        _, protein_name = os.path.split(active_file)
        self.protein_name, _ = os.path.splitext(protein_name)
        self.uniprot_id = protein_to_uniprot.get(self.protein_name, self.protein_name)
        self.xgb_model_file = os.path.join(self.model_dir, f"{self.uniprot_id}.xgb")
        self.rf_model_file = os.path.join(self.model_dir, f"{self.uniprot_id}.rf")
        self.nn_model_file = os.path.join(self.model_dir, f"{self.uniprot_id}.mlp")
        self.svc_model_file = os.path.join(self.model_dir, f"{self.uniprot_id}.svc")
        self.result_file = os.path.join(self.output_dir, f"{self.uniprot_id}_results.json")
        self.results = dict()
        self.refit = 'f1_score'
        self.scoring = {'auc_score': 'roc_auc',
                        'precision_score': make_scorer(metrics.precision_score),
                        'recall_score': make_scorer(metrics.recall_score),
                        'accuracy_score': make_scorer(metrics.accuracy_score)
                        }
        self.fold = 5
        
    def get_data(self):
        try:
            self.actives = pd.read_csv(self.active_file, header=None)
            self.decoys = pd.read_csv(self.decoy_file, nrows=10*len(self.actives), header=None)
        except Exception as e:
            print(str(e))
            return False
                
        actives_x = self.actives.iloc[:, 1:].values
        actives_y = np.ones(len(actives_x))
        decoys_x = self.decoys.iloc[:, :].values
        decoys_y = np.zeros(len(decoys_x))
        self.x = np.concatenate((actives_x, decoys_x)).astype(np.float16)
        self.y = np.concatenate((actives_y, decoys_y)).astype(np.float16)

        # Split the data into train and test
        self.x_train, self.x_test, self.y_train, self.y_test = train_test_split(self.x, self.y, stratify=self.y, test_size=self.test_split, random_state=self.random_state)
        return True
    
    def get_report(self, clf):
        # Save the report
        y_pre = clf.predict(self.x_test)
        y_pro = clf.predict_proba(self.x_test)[:, 1]
        results = dict()
        results['protein_name'] = self.protein_name
        results['uniprot_id'] = self.uniprot_id
        results['roc_auc'] = metrics.roc_auc_score(self.y_test, y_pro)
        results['accuracy'] = metrics.accuracy_score(self.y_test, y_pre)
        results['f1_score'] = metrics.f1_score(self.y_test, y_pre, average='weighted')
        results['cohen_kappa'] = metrics.cohen_kappa_score(self.y_test, y_pre)
        results['mcc'] = metrics.matthews_corrcoef(self.y_test, y_pre)
        results['data_info'] = {"train_count": len(self.x_train), 
                                "test_count": len(self.x_test),
                                "actives_count": len(self.actives), 
                                "decoys_count": len(self.decoys)}
        
        return results

    def write_results(self):
        # Save the report
        with open(self.result_file, 'w') as f:
            json.dump(self.results, f)
    
    def train_xgb(self):
        ratio = float(len(self.actives)) / len(self.decoys)
        clf = xgb.XGBClassifier(max_depth=7,
                               n_jobs=mp.cpu_count(),
                               min_child_weight=1,
                               learning_rate=0.5,
                               n_estimators=1000,
                               silent=True,
                               objective='binary:logistic',
                               gamma=0,
                               max_delta_step=0,
                               subsample=1,
                               colsample_bytree=1,
                               colsample_bylevel=1,
                               reg_alpha=0,
                               reg_lambda=0,
                               scale_pos_weight=ratio,
                               seed=1,
                               missing=None)

        clf.fit(self.x_train, self.y_train, 
                eval_metric=['error', 'logloss'], 
                verbose=False,
                eval_set=[(self.x_train, self.y_train), (self.x_test, self.y_test)], 
                early_stopping_rounds=20)
        
        self.xgb_model = clf
        self.results['xgb'] = self.get_report(clf)    
        # Save the model
        joblib.dump(clf, self.xgb_model_file)
    
    def train_nn(self):
        # Neural Network
        def get_hidden_layers():
            import itertools
            x = [64, 128, 256]
            hl = []

            for i in range(1, len(x)):
                hl.extend([p for p in itertools.product(x, repeat=i+1)])

            return hl

        classifier_nn = MLPClassifier(solver='adam', alpha=1e-5, early_stopping=True, random_state=self.random_state, verbose=classifier_loglevel)
        hidden_layer_sizes = get_hidden_layers()
        parameters_nn = {'hidden_layer_sizes': hidden_layer_sizes}
        gridsearch_nn = GridSearchCV(classifier_nn, parameters_nn, pre_dispatch='n_jobs', scoring=self.scoring, cv=self.fold, refit=self.refit, n_jobs=-1, verbose=gridsearch_loglevel)
        gridsearch_nn.fit(self.x_train, self.y_train)
        self.nn_model = gridsearch_nn.best_estimator_
        self.results['mlp'] = self.get_report(gridsearch_nn.best_estimator_)
        joblib.dump(gridsearch_nn.best_estimator_, self.nn_model_file)
        
    def train_svc(self):
        # Support Vector Machine
        classifier_sv = SVC(class_weight='balanced', kernel='linear', probability=True, random_state=self.random_state, verbose=classifier_loglevel)
        parameters_sv = {'C': [0.1, 1.0, 10, 100, 1000], 'gamma':[0.1, 1, 10, 100, 1000, 'auto']}
        gridsearch_sv = GridSearchCV(classifier_sv, parameters_sv, pre_dispatch='n_jobs', scoring=self.scoring, cv=self.fold, refit=self.refit, n_jobs=-1, verbose=gridsearch_loglevel)
        gridsearch_sv.fit(self.x_train, self.y_train)
        self.svc_model = gridsearch_sv.best_estimator_
        self.results['svc'] = self.get_report(gridsearch_sv.best_estimator_)
        joblib.dump(gridsearch_sv.best_estimator_, self.svc_model_file)

    def train_rf(self):
        classifier_rf = RandomForestClassifier(class_weight='balanced',random_state=self.random_state, verbose=classifier_loglevel)
        parameters_rf = {'n_estimators':[i for i in range(100, 1000, 50)]}
        gridsearch_rf = GridSearchCV(classifier_rf, parameters_rf, pre_dispatch='n_jobs', scoring=self.scoring, cv=self.fold, refit=self.refit, n_jobs=-1, verbose=gridsearch_loglevel)
        gridsearch_rf.fit(self.x_train, self.y_train)
        self.rf_model = gridsearch_rf.best_estimator_
        self.results['rf'] = self.get_report(gridsearch_rf.best_estimator_)
        joblib.dump(gridsearch_rf.best_estimator_, self.rf_model_file)
        
    
    def train_models(self):
        # get data
        if not self.get_data(): 
            print('ERROR: DATA GENERATION FAILED!!')
            return 
        
        # Discard if samples are not enough
        if len(self.actives) < self.sample_threshold: print(f"ERROR: NOT ENOUGH ({len(self.actives)} < {self.sample_threshold}) SAMPLES FOUND!"); return

        print('Training xgboost..')    
        self.train_xgb()
        print('Training random forest model..')    
        self.train_rf()
        print('Training support vector classifier..')    
        self.train_svc()
        print('Training neural network..')    
        self.train_nn()
        
        self.write_results()


if __name__=="__main__":
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(description="Training Ligandnet")
    parser.add_argument('--actives', action='store', dest='actives', required=True,
                        help='actives file (csv)')
    args = vars(parser.parse_args())
    
    actives_file = args['actives']
    _, protein = os.path.split(actives_file)
    protein, _ = os.path.splitext(protein)
    uniprot_id = protein_to_uniprot.get(protein, protein)
    if os.path.isfile(f'reports/{uniprot_id}_results.json'):
        print("ERROR: REPORTS FILE FOUND. MODELS ARE ALREADY TRAINED.")
        sys.exit()
    else:
        print(f'Training {uniprot_id} model')
        
    decoys_file = os.path.join(decoy_dir, f"decoys_{protein}.csv")
    if not os.path.isfile(decoys_file):
        print(f"ERROR: DECOY FILE {decoys_file} NOT FOUND")
        sys.exit()
    t = Train(actives_file, decoys_file)
    t.train_models()
