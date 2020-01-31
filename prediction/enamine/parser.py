# Script for parsing the outputs of Enamine 2M compounds evaluation
# Md Mahmudulla Hassan
# Last modified: 10/17/2018

import pandas as pd
from glob import glob
import numpy as np

advanced_files = glob("output/*/Enamine_advance*", recursive=True)
hts_files = glob("output/*/Enamine_hts*", recursive=True)


def process(files, name):
    output_df =  pd.DataFrame()

    for _file in files:
        df = pd.read_csv(_file)
        df = df.loc[df.prediction > 0.99]
        df = df.loc[df.prediction < 1.0] # Ignore SVMs
        output_df = pd.concat([output_df, df])
    labels = ['protein', 'id', 'smiles', 'prediction']
    output_df.to_csv(name + "_hits.csv", index=None, columns=labels)
    output_df.groupby("protein")['id'].count().reset_index().to_csv(name + "_counts.csv", header=["protein", "count"], index=None)

process(advanced_files, "enamine_advanced")
process(hts_files, "enamine_hts")

