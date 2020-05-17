# Script for predicting ligand activity using LigandNet models
# Author: Md Mahmudulla Hassan
# Department of Computer Science and School of Pharmacy, UTEP
# Last modified: 02/20/2020

from . import LigandNet
import time
import os
import argparse
import errno


if __name__ == "__main__":
    start = time.time()
    parser = argparse.ArgumentParser(
        description="Ligand activity prediction using LigandNet"
    )
    parser.add_argument("--sdf", action="store", dest="sdf", help="SDF file location")
    parser.add_argument(
        "--smiles", action="store", type=str, dest="smiles", help="SMILES"
    )
    #     parser.add_argument('--out', action='store', dest='out',
    #                         required=False, help='Output directory')
    parser.add_argument(
        "--confidence",
        action="store",
        dest="confidence",
        type=float,
        default=0.50,
        help="Minimum confidence to consider for prediction. Default is 0.5",
    )

    args = parser.parse_args()

    if not (args.smiles or args.sdf):
        parser.error("No input found. Provide --smiles or --sdf")

    print("Loading the LigandNet models ...")
    l = LigandNet()

    if args.sdf is not None:
        if not os.path.isfile(args.sdf):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.sdf)

        results = l.get_prediction(args.sdf, "sdf", args.confidence)
        print(results)

    if args.smiles is not None:
        results = l.predict(args.smiles, "smiles", args.confidence)
        print(results)
