#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import

# VERSION NUMBER
__version__ = "1.0.0"
SUPPORTED_EXTENSIONS = set(('.smi','.txt'))

from glob import glob
from argparse import ArgumentParser
import numpy as np
import scipy as sp
import pandas as pd
import os, sklearn, sys
from openbabel import openbabel

# rdkit
import rdkit
from rdkit import Chem
import rdkit.Chem.GraphDescriptors as graph

# external complexity metrics
# SA_Score (Ertl 2009)
try: from .metrics import sa_score
except: from metrics import sa_score

def main():

    # Get command line inputs. Use -h to list all possible arguments and default values
    parser = ArgumentParser()
    parser.add_argument("-addH", dest="ADDH", action="store_true", default=False,
                        help="Add H atoms to all structures before computing complexity score")

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

    for file in files:
        print('o   Reading', file)
        for 

if __name__ == "__main__":
    main()

