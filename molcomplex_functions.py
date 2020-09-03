#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	    	  used for genrating a scores     	    #
#####################################################.

import numpy as np
import pandas as pd
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

# Balaban J Score (Chem. Phys. Lett. 1982, 89, 399-404
def get_balaban_score(mols):
    balaban_scores = []
    for i, mol in enumerate(mols):
        try:
            score = graph.BalabanJ(mol)
        except:
             score = np.nan

        balaban_scores.append(score)
    return balaban_scores

# Kier's alpha-modified shape indices
def get_hallkieralpha_score(mols):
    hkalpha_scores = []
    for i, mol in enumerate(mols):
        try:
            score = graph.HallKierAlpha(mol)
        except:
             score = np.nan

        hkalpha_scores.append(score)
    return hkalpha_scores

#  Bonchev & Trinajstic's information content of the coefficients of the characteristic
# polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule (J. Chem. Phys. 1977, 67, 4517-4533)
def get_ipc_score(mols):
    IPC_scores = []
    for i, mol in enumerate(mols):
        try:
            score = graph.Ipc(mol)
        except:
             score = np.nan

        IPC_scores.append(score)
    return IPC_scores

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
