# coding: utf-8

import json
import os
import pandas as pd

results = json.load(open("results.json", 'r'))
plist = json.load(open("plist.json", "r"))
proteins = [os.path.splitext(i)[0] for i in list(plist.values())[0]]
output = {}
for k, v in results.items():
    for p in proteins:
        if p in v.keys():
            if p in output.keys():
                output[p][k] = v[p]
            else:
                output[p] = {k: v[p]}

with open('results_plist.json', 'w') as f:
    json.dump(output, f)


pd.DataFrame.from_dict(output).to_csv("results_plist.csv")

with open("results_plist.txt", 'w') as f:
    for k, v in output.items():
        f.write(f"{k}\n")
        for p, q in v.items():
            f.write(f"\t{p}: {q}\n")

