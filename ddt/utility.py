# Utility script for feature generation
# Md Mahmudulla Hassan
# The University of Texas at El Paso
# Last Modified: 05/13/2019
from __future__ import print_function, absolute_import
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
import shutil
import subprocess
import errno
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions
import base64
import numpy as np
#import cairosvg

current_dir = os.path.dirname(os.path.realpath(__file__))
MAYACHEMTOOLS_DIR = os.path.join(current_dir, "mayachemtools")

def smiles_to_sdf(smiles):
    # Try to get the rdkit mol
    mol = Chem.MolFromSmiles(smiles)
    # Compute 2D coordinates
    AllChem.Compute2DCoords(mol)
    mol.SetProp("smiles", smiles)
    temp_dir = tempfile.mkdtemp()
    sdf_filepath = os.path.join(temp_dir, "temp.sdf")
    w = Chem.SDWriter(sdf_filepath)
    w.write(mol)
    w.flush()
    return sdf_filepath

class SmilesToImage:
    def __init__(self, smiles):
        self.smiles = smiles
        self.temp_dir = tempfile.mkdtemp()
        self.png_file = os.path.join(self.temp_dir, "mol.png")
        self.svg_file = os.path.join(self.temp_dir, "mol.svg")

    def toPNG(self):
        # Set the drawing options
        DrawingOptions.atomLabelFontSize = 55
        DrawingOptions.dotsPerAngstrom = 100
        DrawingOptions.bondLineWidth = 3.0
        
        # Conver the SMILES into a mol object
        m = Chem.MolFromSmiles(self.smiles)
        # Calculate the coordinates
        AllChem.Compute2DCoords(m)
        # Draw the mol
        Draw.MolToFile(m, self.svg_file)
        # Convert the svg to png (for high quality image)
        #cairosvg.svg2png(url=self.svg_file, write_to=self.png_file)
        # Convert into binary and return
        binary_image = None
        with open(self.png_file, "rb") as f:
            binary_image = base64.b64encode(f.read())
            shutil.rmtree(self.temp_dir)

        return binary_image


class FeatureGenerator:
    
    def __init__(self):
        self.sdf_filepath = None
        self.smiles = None
        
    def load_smiles(self, smiles):
        self.smiles = smiles
        self.sdf_filepath = smiles_to_sdf(self.smiles)
    
    def load_sdf(self, sdf_filepath):
        self.sdf_filepath = sdf_filepath    
    
    def extract_tpatf(self, save_csv=False):
        features = []
        script_path = os.path.join(MAYACHEMTOOLS_DIR, "bin/TopologicalPharmacophoreAtomTripletsFingerprints.pl")
        # Generate the TPATF features
        # Check if the sdf file exists
        if not os.path.isfile(self.sdf_filepath):
            print("SDF file not found")
            return None
        
        # Extract tpatf features
        temp_dir = tempfile.mkdtemp()    
        temp_file = os.path.join(temp_dir, "temp")
        # Use below command to use the compound name as compound id (selleck compounds sdf has it)
        #command = "perl " + script_path + " -r " + temp_file + " --DataFieldsMode Specify --DataFields Name --AtomTripletsSetSizeToUse FixedSize -v ValuesString -o " + self.sdf_filepath
        # Use below command to use the available compound id in the sdf file (most sdf files have it)
        command = "perl " + script_path + " -r " + temp_file + " --CompoundIDMode MolnameOrLabelPrefix --AtomTripletsSetSizeToUse FixedSize -v ValuesString -o " + self.sdf_filepath
        #command = "perl " + script_path + " -r " + temp_file + " --DataFieldsMode CompoundID --CompoundIDMode MolnameOrLabelPrefix --CompoundID Cmpd --CompoundIDLabel MolID --AtomTripletsSetSizeToUse FixedSize -v ValuesString -o " + self.sdf_filepath
        os.system(command)
        output_csv = temp_file + ".csv"
        if not os.path.isfile(output_csv):
            print("ERROR: TPATF features wasn't extracted")
            return None

        if save_csv: shutil.copy(output_csv, os.getcwd())
        
        compound_list = []
        content = []
        
        #TODO: Add warning with line number that fails to load because of invalid decode error
        with open(output_csv, 'r') as f:
            line = f.readline()
            content.append(line)
            while line:
                try:
                    line = f.readline()
                    if line: content.append(line) # Checks emptry string
                except:
                    continue
                    
        content = [c.replace('"', '') for c in content] # remove (") from the content
        content = [c.split(';') for c in content] # split features from the other information
        # Separate the compound features
        for con in content[1:]: # First item doesn't have features
            if len(con) < 6:
                print("Inconsistent feature entry. Ignoring {}".format(con[0].split(',')[0]))
                continue
            compound_list.append(con[0].split(',')[0])
            features.append([int(i) for i in con[5].split(" ")])

        #with open(output_csv, 'r') as f:
        #    for line in f.readlines():
        #        #if "Cmpd" in line:
        #        if line[1:5] == "Cmpd":
        #            #compound_list.append(list[1:
        #            line = line.split(';')[5].replace('"','')
        #            features.append([int(i) for i in line.split(" ")])

        # Clean up the temporary files
        shutil.rmtree(temp_dir)
        return compound_list, np.array(features).reshape((-1, 2692)).astype(np.float32)


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="An utility script for small tasks.")
    parser.add_argument('--smiles', action='store', dest='smiles', required=False, help='SMILES string as "SMILES"')
    parser.add_argument('--sdf', action='store', dest='sdf', required=False, help='SDF file location')
    parser.add_argument('--to', action='store', dest='to', required=False, help='Save features to .npy')
    parser.add_argument('--csv', action='store_true', dest='csv', required=False, help='Get the features in a csv file')
    parse_dict = vars(parser.parse_args())
    
    # Example: Extracting TPATF features
    ft = FeatureGenerator()
    if parse_dict['smiles'] is not None:
        ft.load_smiles(parse_dict['smiles'])
    elif parse_dict['sdf'] is not None:
        ft.load_sdf(parse_dict['sdf'])
        compounds, features = ft.extract_tpatf(True) if parse_dict['csv'] else ft.extract_tpatf()
        if parse_dict['to'] is not None:
            np.save(parse_dict['to'], features)
            with open(parse_dict['to'] + '.txt', 'w') as f:
              f.writelines(','.join([c for c in compounds]))
        else:
            print(compounds, len(features))
