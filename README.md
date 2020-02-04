# LigandNet
LigandNet, a tool which combines different machine learning models into one platform for the prediction of the state of the ligands either actives or inactives for a particular proteins.

# Setup
Create a conda environment using `environment.yml`. Run the following
```bash
conda env create -f environment.yml
```

# Run
Use `ligandnet.py` to run predictions. To see the available options, run `python ligandnet.py --help` which shows the following:

```bash
usage: ligandnet.py [-h] [--sdf SDF] [--smiles SMILES]
                    [--confidence CONFIDENCE]

Ligand activity prediction using LigandNet

optional arguments:
  -h, --help            show this help message and exit
  --sdf SDF             SDF file location
  --smiles SMILES       SMILES
  --confidence CONFIDENCE
                        Minimum confidence to consider for prediction. Default
                        is 0.5
```

For example, `python ligandnet.py --smiles CCCC` will run all the LigandNet models on the compound `CCCC`. For an sdf file as input, run `python ligandnet.py --sdf samples/AAAAML.xaa.sdf`. The parameter `confidence` is the minimum probability for which a model will consider a ligand as an active.

# Decoys
To get the decoys used for training the LigandNet models, run 

```bash
1. bash get_decoys.sh
2. tar xvf decoys.tar.gz
```

# Web server
An web interface for ligand activity prediction using the LigandNet models is available at [LigandNet](https://drugdiscovery.utep.edu/ligandnet)
