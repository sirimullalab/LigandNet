# Script for fetching the gene names from Pharos dataset for corresponding targets
# Md Mahmudulla Hassan
# Last modified: 10/19/2018

import json
from urllib.request import urlopen
from urllib.parse import quote
import sys
import pandas as pd
pharos = 'https://pharos.nih.gov/idg/api/v1/targets/'
#curl -X GET "https://pharos.nih.gov/idg/api/v1/targets/search?q=sodium%20dependent%20dopamine%20transporter&top=10"

def fetch_ligand(target):
    req = urlopen(pharos + target)
    target = json.loads(req.read())
    # print(json.dumps(target, indent=4, separators=(',', ': ')))
    # retrieve ligand links for this target
    req = urlopen(pharos+'{0}'.format(target['id'])
                  +'/links(kind=ix.idg.models.Ligand)')
    link = json.loads(req.read())
    # print(json.dumps(link, indent=4, separators=(',', ': ')))
    
    for l in link:
        name = ""
        for p in l['properties']:
            if p['label'] == 'IDG Ligand' or p['label'] == 'IDG Drug':
                name = p['term']
                break
        
        req = urlopen(l['href']+'/properties(label=CHEMBL Canonical SMILES)')
        ligand = json.loads(req.read())
        
        print(ligand[0]['text'] + '\t' + name, end=' ')
        # do another pass through the properties to extract ligand activity
        for p in l['properties']:
            if p['label'] == 'Ki' or p['label'] == 'Kd' or p['label'] == 'IC50' or p['label'] == 'EC50' or p['label'] == 'AC50' or p['label'] == 'Potency':
                print('\tp' + p['label']+'={0}'.format(p['numval']), end=' ')
        print()
    

def fetch_gene(target, query="gene"):
    # replace '_' with space
    target = target.replace("_", " ")
    print("Fetching for " + target, end=' ')
    q = pharos + "search?q=" + quote(target)
    req = urlopen(q)
    data = json.loads(req.read().decode())
    try:
        content = data['content'][0]
        gene = content[query]
        print(" -> ", gene)
        return gene
    except Exception as e:
        print("ERROR in " + target)
        print(data)
        print(str(e))
        return " "

def fetch_gene_driver(_file):
    df = pd.read_csv(_file)
    df['gene'] = df.protein.apply(fetch_gene)
    # print(df[['protein', 'gene']])
    print("All genes are fetched and also written into a csv file.")
    df.to_csv(_file + "_genes.csv", index=None)


if __name__ == "__main__":
    #fetch_gene_driver("../enamine_advanced_counts.csv")
    #fetch_gene_driver("../enamine_hts_counts.csv")
    #fetch_gene_driver("/data2/mhassan/LigandNet/models/all_models_types.csv")
    fetch_gene_driver("../chembridge/counts.csv")
