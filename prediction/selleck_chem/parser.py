import json

results = json.load(open("py27_results.json", "r"))
results.update(json.load(open("py35_results.json", "r")))

with open("active_compounds.txt", "w") as f:
    for k, v in results.items():
        f.write("PROTEIN: {}\nACTIVE COMPOUNDS: {}\n\n".format(k, v))
        # break
