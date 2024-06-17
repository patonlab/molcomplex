#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import

# VERSION NUMBER
__version__ = "1.1.0"
SUPPORTED_EXTENSIONS = set(('.smi','.txt'))

from argparse import ArgumentParser
from glob import glob

import numpy as np
import pandas as pd
import os, time

from molcomplex.descriptors import *
from molcomplex.complex_object import *
from molcomplex.complex_funcs import *

def main():
    start = time.time()
    # Get command line inputs. Use -h to list all possible arguments and default values
    parser = ArgumentParser()
    parser.add_argument("--addH", dest="addH", action="store_true", default=False,
                        help="Add H atoms to all structures before computing complexity score")
    parser.add_argument("--csv", dest="csv", action="store_true", default=False,
                        help="Export data to .CSV formatted file")
    parser.add_argument("-f","--file", dest="file", default=[],
                        help="Input file for calculations of descriptors",nargs='*')
    parser.add_argument("-o","--output", dest="output", default='output',
                        help="Output filename for calculations of descriptors to be stored in CSV")
    parser.add_argument("--twc", dest="twc", action="store_true", default=False,
                        help="Rucker's total walk count")
    parser.add_argument("--retro", dest="retro", action="store_true", default=False,
                        help="Perform retrosynthesis on the list of molecules provided")
    parser.add_argument("--nbonds", dest="nbonds", action="store", default=1,
                        help="Number of bonds to break in retrosynthesis disconnections")

    # Parse Arguments
    (options, args) = parser.parse_known_args()

    # File(s) to be parsed
    files = []
    for elem in options.file:
        if os.path.splitext(elem)[1].lower() in SUPPORTED_EXTENSIONS:  # Look for file names
            for file in glob(elem): files.append(file)
        else:
            print(' ! ', elem, 'is not currently supported !')

    # this assumes that we are dealing with lists of smiles
    mol_smiles = []

    for file in files:
        print('o   Reading', file)
        smifile = open(file)

        for line in smifile:
            toks = line.split()
            if len(toks) != 0:
                smi = toks[0]
                mol_smiles.append(smi)
    
    if options.retro:
        df = parse_contents(mol_smiles, options.nbonds)
        mol_smiles = df.SMILES

    MOL_DATA = mol_complex(mol_smiles, options.twc)
    if options.retro:
        MOL_DATA = MOL_DATA.merge(df, on=['SMILES'])

    total_time = time.time() - start

    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(MOL_DATA)

    print("Num Descriptors: ",len(MOL_DATA.columns) - 1)
    print("Total Time: ",total_time)

    if options.csv:
        MOL_DATA.to_csv(f"molcomplex_{options.output}.csv",index=False)

if __name__ == "__main__":
    main()
