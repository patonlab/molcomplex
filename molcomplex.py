#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import

# VERSION NUMBER
__version__ = "1.0.0"
SUPPORTED_EXTENSIONS = set(('.smi','.txt'))

from argparse import ArgumentParser
from glob import glob

import numpy as np
import pandas as pd
import os, sys, time
from concurrent import futures
max_workers = 2

from molcomplex_functions import *
from complex_object import *

def main():
    start = time.time()
    # Get command line inputs. Use -h to list all possible arguments and default values
    parser = ArgumentParser()
    parser.add_argument("--addH", dest="addH", action="store_true", default=False,
                        help="Add H atoms to all structures before computing complexity score")
    parser.add_argument("--csv", dest="csv", action="store_true", default=False,
                        help="Export data to .CSV formatted file")
    # Parse Arguments
    (options, args) = parser.parse_known_args()

    # File(s) to be parsed
    files = []
    args = sys.argv[1:]
    for elem in args:
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
                pieces = smi.split('.')
                if len(pieces) > 1:
                    smi = max(pieces, key=len) #take largest component by length
                    print("Taking largest component: %s\t%s" % (smi,name))

                # append all smiles to the mol_smiles list
                mol_smiles.append(smi)
    
    MOL_DATA = mol_complex(mol_smiles)

    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(MOL_DATA)
    total_time = time.time() - start
    
    print("Num Descriptors: ",len(MOL_DATA.columns) - 1)
    print("Total Time: ",total_time)
    if options.csv:
        MOL_DATA.to_csv("molcomplex.csv",index=False)
    
if __name__ == "__main__":
    main()
