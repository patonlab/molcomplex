#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	        used for genrating a descriptors     	#
#####################################################.

import os
import numpy as np

from rdkit import Chem
from rdkit.Chem.AtomPairs import Pairs,Torsions
import rdkit.Chem.GraphDescriptors as graph
from openbabel import openbabel
from mordred import CPSA, KappaShapeIndex, McGowanVolume, MoeType, VdwVolumeABC, ZagrebIndex

# external complexity metrics
try:
	from .metrics import sa_score
	from .metrics import boettcher
	from .metrics import standalone_model_numpy
	from .metrics import proudfoot
	from .metrics import rucker_twc
except:
	from metrics import sa_score
	from metrics import boettcher
	from metrics import standalone_model_numpy
	from metrics import proudfoot
	from metrics import rucker_twc


# Bertz Complexity Score (JACS 1981, 103, 3241-3243)
def get_bertz_score(mol):
	try:
		score = graph.BertzCT(mol)
	except:
		score = np.nan
	return score

# Balaban J Score (Chem. Phys. Lett. 1982, 89, 399-404
def get_balaban_score(mol):
	try:
		score = graph.BalabanJ(mol)
	except:
		score = np.nan
	return score

# Kier's alpha-modified shape indices
def get_hallkieralpha_score(mol):
	try:
		score = graph.HallKierAlpha(mol)
	except:
		score = np.nan
	return score

#  Bonchev & Trinajstic's information content of the coefficients of the characteristic
# polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule (J. Chem. Phys. 1977, 67, 4517-4533)
def get_ipc_score(mol):
	try:
		score = graph.Ipc(mol)
	except:
		score = np.nan
	return score

# Ertl SA_Score (J. Cheminform. 2009, 1, 8)
def get_sa_score(mol):
	# read fragment scores from file
	sa_score.readFragmentScores("fpscores")
	
	try:
		score = sa_score.calculateScore(mol)
	except:
		score = np.nan
	return score

# Boettcher Score (J. Chem. Inf. Model. 2016, 56, 3, 462–470)
def get_boettcher_score(mol):
	#needs OBMol input to compute automorphism with openbabel
	smi = Chem.MolToSmiles(mol)
	obConversion = openbabel.OBConversion()
	obConversion.SetInAndOutFormats("smi", "smi")
	bottch = boettcher.BottchScore("False")

	try:
		obmol = openbabel.OBMol()
		obConversion.ReadString(obmol, smi)
		score=bottch.score(obmol)
	except:
		score = np.nan

	return score

#  SCScore (J. Chem. Inf. Model. 2018, 58, 2, 252)
def get_scscore(mol,model):
	smi = Chem.MolToSmiles(mol)
	try:
		(smi, score) = model.get_score_from_smi(smi)
	except:
		score = np.nan

	return score

def get_rucker_twc(mol):
	try:
		score = rucker_twc.twc(mol)
	except:
		score = np.nan

	return score

def get_proudfoot_index(mol):
	try:
		score = proudfoot.proudfoot_index(mol)[0]
	except:
		score = np.nan

	return score

''' Taken from Merck Paper descriptors '''

''' 1. Descriptor complex - AP and TT '''
def DESCRIPTORCOMPLEXITY_UNIQUEAP(mol):
	try:
		mol_noHs = Chem.RemoveHs(mol)
		mol_ap_fp = Pairs.GetAtomPairFingerprint(mol_noHs)
		num_uniq_ap = len(mol_ap_fp.GetNonzeroElements())
	except:
		num_uniq_ap = np.nan

	return num_uniq_ap

def DESCRIPTORCOMPLEXITY_UNIQUETT(mol):
	# need to check as the example from merck as 19 unique but actual unique is 13.
	
	try:
		mol_noHs = Chem.RemoveHs(mol)
		mol_tt_fp = Torsions.GetTopologicalTorsionFingerprint(mol_noHs)
		num_uniq_tt = len(mol_tt_fp.GetNonzeroElements())
	except:
		num_uniq_tt = np.nan

	return num_uniq_tt

def DESCRIPTORCOMPLEXITY_TOTALAP(mol):	
	try:
		mol_noHs = Chem.RemoveHs(mol)
		num_noHs = mol_noHs.GetNumAtoms()
		num_tot_AP = (num_noHs*(num_noHs - 1))/2
	except:
		num_tot_AP= np.nan

	return num_tot_AP

def DESCRIPTORCOMPLEXITY_TOTALTT(mol):
	try:
		num_tot_TT = len(Chem.rdmolops.FindAllPathsOfLengthN(mol,3))
	except:
		num_tot_TT = np.nan

	return num_tot_TT


def DESCRIPTORCOMPLEXITY_APCOMPLEX(mol):
	num_uniq = DESCRIPTORCOMPLEXITY_UNIQUEAP(mol)
	num_tot = DESCRIPTORCOMPLEXITY_TOTALAP(mol)

	if num_tot == 0: apc = 0
	else: apc = num_uniq/num_tot

	return apc


def DESCRIPTORCOMPLEXITY_TTCOMPLEX(mol):
	num_uniq = DESCRIPTORCOMPLEXITY_UNIQUETT(mol)
	num_tot = DESCRIPTORCOMPLEXITY_TOTALTT(mol)
	
	if num_tot == 0: ttc = 0
	else: ttc = num_uniq/num_tot

	return ttc

''' 2. MOE_2D '''


''' 3. SP3CARBONS '''

def SP3CARBONS_TOTALATOM_COUNT(mol):
	try:
		totatoms = mol.GetNumAtoms()
	except:
		totatoms = np.nan

	return totatoms


def SP3CARBONS_TOTALCARBON_COUNT(mol):
	try:
		totcar =  0
		for atom in mol.GetAtoms():
			if atom.GetSymbol() == 'C':
				totcar += 1
	except:
		totcar = np.nan

	return totcar


def SP3CARBONS_CAR_ALLATOM_RATIO(mol):
	numcar = SP3CARBONS_TOTALCARBON_COUNT(mol)
	numatoms  = SP3CARBONS_TOTALATOM_COUNT(mol)

	if numatoms == 0: cratio = 0
	else: cratio = numcar / numatoms

	return cratio

######################

def SP3CARBONS_CHIRAL_ALLATOM_RATIO(mol):
	numchiral = SP3CARBONS_CHIRAL_COUNT(mol)
	numatoms  = SP3CARBONS_TOTALATOM_COUNT(mol)

	if numatoms == 0: chiratio = 0
	else: chiratio = numchiral / numatoms	
	
	return chiratio


def SP3CARBONS_CHIRAL_ALLCARBON_RATIO(mol):
	numchiral = SP3CARBONS_CHIRAL_COUNT(mol)
	numcarbons  = SP3CARBONS_TOTALCARBON_COUNT(mol)

	if numcarbons == 0: chiratio = 0
	else: chiratio = numchiral / numcarbons

	return chiratio


def SP3CARBONS_CHIRAL_COUNT(mol):
	try:
		c_chiral = 0
		chiralcenters = Chem.FindMolChiralCenters(mol,includeUnassigned=True)
		for list in chiralcenters:
			if mol.GetAtomWithIdx(list[0]).GetSymbol() == 'C':
				c_chiral += 1
	except:
		c_chiral = np.nan

	return c_chiral

#########################

def SP3CARBONS_CSP2_COUNT(mol):
	try:
		csp2 =  0
		for atom in mol.GetAtoms():
			if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP2:
				csp2 += 1
	except:
		csp2 = np.nan
	return csp2

def SP3CARBONS_CSP2_ALLATOM_RATIO(mol):
	numcsp2 = SP3CARBONS_CSP2_COUNT(mol)
	numatoms  = SP3CARBONS_TOTALATOM_COUNT(mol)

	if numatoms == 0: csp2ratio = 0
	else: csp2ratio = numcsp2 / numatoms

	return csp2ratio


def SP3CARBONS_CSP2_ALLCARBON_RATIO(mol):
	numcsp2 = SP3CARBONS_CSP2_COUNT(mol)
	numcar  = SP3CARBONS_TOTALCARBON_COUNT(mol)

	if numcar == 0: csp2ratio = 0
	else: csp2ratio = numcsp2 / numcar

	return csp2ratio

##########################

def SP3CARBONS_CSP3_COUNT(mol):
	try:
		csp3 =  0
		for atom in mol.GetAtoms():
			if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP3:
				csp3 += 1
	except:
		csp3 = np.nan

	return csp3


def SP3CARBONS_CSP3_ALLATOM_RATIO(mol):
	numcsp3 = SP3CARBONS_CSP3_COUNT(mol)
	numatoms  = SP3CARBONS_TOTALATOM_COUNT(mol)

	if numatoms == 0: csp3ratio = 0
	else: csp3ratio = numcsp3 / numatoms

	return csp3ratio


def SP3CARBONS_CSP3_ALLCARBON_RATIO(mol):
	numcsp3 = SP3CARBONS_CSP3_COUNT(mol)
	numcar  = SP3CARBONS_TOTALCARBON_COUNT(mol)

	if numcar == 0: csp3ratio = 0
	else: csp3ratio = numcsp3 / numcar

	return csp3ratio

###########################

def SP3CARBONS_CSP_COUNT(mol):
	try:
		csp =  0
		for atom in mol.GetAtoms():
			if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP:
				csp += 1
	except:
		csp = np.nan
	return csp


def SP3CARBONS_CSP_ALLATOM_RATIO(mol):
	numcsp = SP3CARBONS_CSP_COUNT(mol)
	numatoms  = SP3CARBONS_TOTALATOM_COUNT(mol)

	if numatoms == 0: cspratio = 0
	else: cspratio = numcsp / numatoms
	
	return cspratio


def SP3CARBONS_CSP_ALLCARBON_RATIO(mol):
	numcsp = SP3CARBONS_CSP_COUNT(mol)
	numcar  = SP3CARBONS_TOTALCARBON_COUNT(mol)
	
	if numcar == 0: cspratio = 0
	else: cspratio = numcsp / numcar

	return cspratio

''' general rdkit descriptors  - rings information'''

def RINGINFO_NUM_ALI_CARBOCYCLE(mol):
	try:
		aliccycle = Chem.rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)
	except:
		aliccycle = np.nan
	return aliccycle


def RINGINFO_NUM_ALI_HETEROCYCLE(mol):
	try:
		alihetcycle = Chem.rdMolDescriptors.CalcNumAliphaticHeterocycles(mol)
	except:
		alihetcycle = np.nan
	return alihetcycle

def RINGINFO_NUM_ALI_RINGS(mol):
	try:
		alirings = Chem.rdMolDescriptors.CalcNumAliphaticRings(mol)
	except:
		alirings= np.nan
	return alirings


def RINGINFO_NUM_ARO_CARBOCYCLE(mol):
	try:
		aroccycle = Chem.rdMolDescriptors.CalcNumAromaticCarbocycles(mol)
	except:
		aroccycle= np.nan
	return aroccycle


def RINGINFO_NUM_ARO_HETEROCYCLE(mol):
	try:
		arohetcycle = Chem.rdMolDescriptors.CalcNumAromaticHeterocycles(mol)
	except:
		arohetcycle= np.nan

	return arohetcycle


def RINGINFO_NUM_ARO_RINGS(mol):
	try:
		arorings = Chem.rdMolDescriptors.CalcNumAromaticRings(mol)
	except:
		arorings= np.nan

	return arorings


def RINGINFO_NUM_BRIDGE_ATOMS(mol):
	try:
		bridge = Chem.rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
	except:
		bridge = np.nan

	return bridge


def RINGINFO_NUM_SPIRO_ATOMS(mol):
	try:
		spiro = Chem.rdMolDescriptors.CalcNumSpiroAtoms(mol)
	except:
		spiro= np.nan

	return spiro

''' RDKit descriptors - present '''

def WIENER_INDEX(mol):
	try:
		wieneridx = 0
		amat = Chem.GetDistanceMatrix(mol)
		num_atoms = mol.GetNumAtoms()
		for i in range(num_atoms):
			for j in range(i+1,num_atoms):
				wieneridx += amat[i][j]
	except:
		wieneridx = np.nan

	return wieneridx


def SMILES_3_2(mol):
	smi = Chem.MolToSmiles(mol)
	try:
		smi32 = len(smi)**(3/2)
	except:
		smi32 = np.nan

	return smi32

''' Taken from pubchem descriptors '''

def PUBCHEM_XLOGP(mol):
	try:
		xlogp = Chem.rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
	except:
		xlogp = np.nan

	return xlogp


def PUBCHEM_TPSA(mol):
	try:
		tpsa = Chem.rdMolDescriptors.CalcTPSA(mol)
	except:
		tpsa = np.nan

	return tpsa


def PUBCHEM_H_BOND_DONOR_COUNT(mol):
	try:
		hbd = Chem.rdMolDescriptors.CalcNumHBD(mol)
	except:
		hbd = np.nan

	return hbd


def PUBCHEM_H_BOND_ACCEPTOR_COUNT(mol):
	try:
		hba = Chem.rdMolDescriptors.CalcNumHBA(mol)
	except:
		hba = np.nan

	return hba


def PUBCHEM_ROTATABLE_BOND_COUNT(mol):
	try:
		rbond = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol,strict=1)
	except:
		rbond = np.nan

	return rbond


def PUBCHEM_HEAVY_ATOM_COUNT(mol):
	try:
		hac = Chem.Lipinski.HeavyAtomCount(mol)
	except:
		hac = np.nan
	
	return hac


def PUBCHEM_ATOM_STEREO_COUNT(mol):
	try:
		asc = Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol)
	except:
		asc = np.nan

	return asc


def PUBCHEM_DEFINED_ATOM_STEREO_COUNT(mol):
	all_count = PUBCHEM_ATOM_STEREO_COUNT(mol)
	undefined_count = PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(mol)

	dasc = all_count-undefined_count

	return dasc


def PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(mol):
	try:
		uasc = Chem.rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)
	except:
		uasc = np.nan
	
	return uasc


def PUBCHEM_BOND_STEREO_COUNT(mol):
	return 0


def PUBCHEM_DEFINED_BOND_STEREO_COUNT(mol):
	return 0


def PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(mol):
	return 0


def PUBCHEM_COVALENT_UNIT_COUNT(mol):
	return 0

''' Mordred provided descriptors '''

def KAPPA_SHAPE_INDEX1(mol):
	"""Defined on page 9 of https://www.epa.gov/sites/production/files/2015-05/documents/moleculardescriptorsguide-v102.pdf"""
	ksi1 = KappaShapeIndex.KappaShapeIndex1()

	try:
		ksi = ksi1(mol)
	except:
		ksi = np.nan

	return ksi


def KAPPA_SHAPE_INDEX2(mol):
	ksi2 = KappaShapeIndex.KappaShapeIndex2()
	try:
		ksi = ksi2(mol)
	except:
		ksi = np.nan
	
	return ksi


def KAPPA_SHAPE_INDEX3(mol):
	ksi3 = KappaShapeIndex.KappaShapeIndex3()
	try:
		ksi = ksi3(mol)
	except:
		ksi = np.nan

	return ksi


def MCGOWAN_VOLUME(mol):
	mv = McGowanVolume.McGowanVolume()
	try:
		mv_num = mv(mol)
	except:
		mv_num = np.nan
	
	return mv_num


def MOE_TYPE_Labute_ASA(mol):
	lasa = MoeType.LabuteASA()
	try:
		lasa_num = lasa(mol)
	except:
		lasa_num = np.nan
	
	return lasa_num


def MOE_TYPE_PEOE_VSA(mol):
	pvsa = MoeType.PEOE_VSA()
	try:
		pvsa_val = pvsa(mol)
	except:
		pvsa_val = np.nan
	
	return pvsa_val


def MOE_TYPE_SMR_VSA(mol):
	svsa = MoeType.SMR_VSA()
	try:
		svsa_val =  svsa(mol)
	except:
		svsa_val= np.nan
	
	return svsa_val


def MOE_TYPE_SLOGP_VSA(mol):
	slpvsa = MoeType.SlogP_VSA()
	try:
		slpvsa_val = slpvsa(mol)
	except:
		slpvsa_val= np.nan
	
	return slpvsa_val


def MOE_TYPE_ESTATE_VSA(mol):
	esvsa = MoeType.EState_VSA()
	try:
		esvsa_val = esvsa(mol)
	except:
		esvsa_val= np.nan
	
	return esvsa_val


def VDW_VOLUME_ABC(mol):
	vvabc = VdwVolumeABC.VdwVolumeABC()
	try:
		vvabc_val = vvabc(mol)
	except:
		vvabc_val= np.nan
	
	return vvabc_val


def ZAGREB_INDEX(mol):
	zi = ZagrebIndex.ZagrebIndex()
	try:
		zi_val = zi(mol)
	except:
		zi_val= np.nan
	
	return zi_val

#template for new descriptors
#
# 	try:
#		score = 
# 	except:
# 		score = np.nan
# 
# return score
