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
try:
    from .metrics import sa_score
    from .metrics import boettcher
    from .metrics import standalone_model_numpy
except:
    from metrics import sa_score
    from metrics import boettcher
    from metrics import standalone_model_numpy

# Bertz Complexity Score (JACS 1981, 103, 3241-3243)
def get_bertz_score(mols):
    bertz_scores = []
    for i, mol in enumerate(mols):
        try:
            score = graph.BertzCT(mol)
        except:
             score = np.nan

        bertz_scores.append(score)
    return bertz_scores

# Ertl SA_Score (J. Cheminform. 2009, 1, 8)
def get_sa_score(mols):
    # read fragment scores from file
    sa_score.readFragmentScores("fpscores")
    SA_Scores = []
    for i, mol in enumerate(mols):
        try:
            score = sa_score.calculateScore(mol)
        except:
            score = np.nan
        SA_Scores.append(score)
    return SA_Scores

# Boettcher Score (J. Chem. Inf. Model. 2016, 56, 3, 462–470)
def get_boettcher_score(mols):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "smi")
    bottch = boettcher.BottchScore("False")
    Boettcher_Scores = []

    for i, smi in enumerate(mols):
        try:
            mol = openbabel.OBMol()
            obConversion.ReadString(mol, smi)
            score=bottch.score(mol)
        except:
            score = np.nan

        Boettcher_Scores.append(score)
    return Boettcher_Scores

#  SCScore (J. Chem. Inf. Model. 2018, 58, 2, 252)
def get_scscore(mols):
    model = standalone_model_numpy.SCScorer()
    model.restore(os.path.join('.', 'models', 'full_reaxys_model_1024bool', 'model.ckpt-10654.as_numpy.json.gz'))

    SC_Scores = []
    for i, smi in enumerate(mols):
        try:
            (smi, score) = model.get_score_from_smi(smi)
        except:
            score = np.nan
        SC_Scores.append(score)
    return SC_Scores


def main():

    # Get command line inputs. Use -h to list all possible arguments and default values
    parser = ArgumentParser()
    parser.add_argument("-addH", dest="addH", action="store_true", default=False,
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

    # this assumes that we are dealing with lists of smiles
    mol_smiles = []

    for file in files:
        print('o   Reading', file)
        smifile = open(file)

        for line in smifile:
            toks = line.split()
            smi = toks[0]

            pieces = smi.split('.')
            if len(pieces) > 1:
                smi = max(pieces, key=len) #take largest component by length
                print("Taking largest component: %s\t%s" % (smi,name))

            # append all smiles to the mol_smiles list
            mol_smiles.append(smi)

    # create a dataframe with a smiles column; the metrics will be added as new columns
    MOL_DATA = pd.DataFrame(mol_smiles, columns=['SMILES'])

    # create a list of rdkit mol objects (required for some metrics)
    mol_objects = [Chem.MolFromSmiles(smi) for smi in mol_smiles]
    if options.addH == True: mol_objects = [Chem.AddHs(mol) for mol in mol_objects]

    # obtain compexity scores for the entire list
    MOL_DATA['SAS'] = get_sa_score(mol_objects)
    MOL_DATA['BERTZ'] = get_bertz_score(mol_objects)
    MOL_DATA['BOETTCHER'] = get_boettcher_score(mol_smiles)
    #MOL_DATA['SCS'] = get_scscore(mol_smiles)


    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(MOL_DATA)

if __name__ == "__main__":
    main()
