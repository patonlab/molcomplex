#/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	    used for genrating a Rucker total  index   #
#####################################################.

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def awc_count(walk,atom):
    if walk == 1:
        return atom.GetDegree()
    else:
        count = 0
        for neigh in atom.GetNeighbors():
            count += awc_count(walk-1, neigh)
        return count

def twc(mol):
    awc = np.zeros((mol.GetNumAtoms(), mol.GetNumAtoms()-1))
    for atom_num, atom in enumerate(mol.GetAtoms()):
        for walk in range(1,mol.GetNumAtoms()):
            awc[atom_num][walk-1] = awc_count(walk,atom)

    mwc_count =  np.sum(awc, axis = 0)
    twc_count = np.sum(mwc_count)/2
    return twc_count
