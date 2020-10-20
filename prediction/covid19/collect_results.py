# coding: utf-8
import pandas as pd
import json
import os
from glob import glob

drugs = pd.read_csv("drugs_all.csv")
result_files = glob("result_files/*.out")
unicode_to_protein = json.load(open('../../uniprot_to_protein_name.json'))
result_dict = {}
for _, row in drugs.iterrows():
    d = row['drug']; s = row['smiles']
    if not os.path.isfile('result_files/' + s + '.out'): continue
    result_dict[d] = {'smiles': s, 'results': json.load(open('result_files/' + s + '.out'))['Cmpd1']}
      
for key, value in result_dict.items():
    output_file = os.path.join('predictions', key + ".csv")
    df = pd.DataFrame.from_dict(value['results'], orient='index', columns=['confidence'])
    df.index.name = 'uniprot_id'; df = df.reset_index()
    df['protein_name'] = df.uniprot_id.apply(lambda x: unicode_to_protein[x])
    df.to_csv(output_file, index=False)
    
