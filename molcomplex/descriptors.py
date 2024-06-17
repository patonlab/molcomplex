#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	        used for genrating a descriptors     	#
#####################################################.

import os
import numpy as np

import networkx as nx
from rdkit import Chem
from rdkit.Chem.AtomPairs import Pairs,Torsions
import rdkit.Chem.GraphDescriptors as graph
from rdkit.Chem import Descriptors
from openbabel import openbabel
openbabel.obErrorLog.SetOutputLevel(0)
from mordred import CPSA, KappaShapeIndex, McGowanVolume, MoeType, VdwVolumeABC, ZagrebIndex
from syba.syba import SybaClassifier
import pkg_resources

#import dbstep.Dbstep as db

# external complexity metrics
try:
	from molcomplex.metrics import sa_score
	from molcomplex.metrics import boettcher
	from molcomplex.metrics import standalone_model_numpy
	from molcomplex.metrics import proudfoot
	from molcomplex.metrics import rucker_twc
	from molcomplex.metrics import sps_score
except:
	from metrics import sa_score
	from metrics import boettcher
	from metrics import standalone_model_numpy
	from metrics import proudfoot
	from metrics import rucker_twc
	from metrics import sps_score


# Bertz Complexity Score (JACS 1981, 103, 3241-3243)
def get_bertz_score(mols):
	bertz_scores = []
	for i, mol in enumerate(mols):
		score = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				score += graph.BertzCT(mol)
			except:
				pass

		bertz_scores.append(score)
	return bertz_scores

# Balaban J Score (Chem. Phys. Lett. 1982, 89, 399-404
def get_balaban_score(mols):
	balaban_scores = []
	for i, mol in enumerate(mols):
		score = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				score += graph.BalabanJ(mol)
			except:
				pass

		balaban_scores.append(score)
	return balaban_scores

# Spatial score based on https://doi.org/10.1021/acs.jmedchem.3c00689.
def get_spatial_score(smis):
	spatial_score = []
	for i, smi in enumerate(smis):
		score = 0
		for s in smi.split('.'):
			try:
				score += sps_score.calculate_score_from_smiles(s)
			except:
				pass

		spatial_score.append(score)
	return spatial_score


# Kier's alpha-modified shape indices
def get_hallkieralpha_score(mols):
	hkalpha_scores = []
	for i, mol in enumerate(mols):
		score = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				score += graph.HallKierAlpha(mol)
			except:
				pass

		hkalpha_scores.append(score)
	return hkalpha_scores

#  Bonchev & Trinajstic's information content of the coefficients of the characteristic
# polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule (J. Chem. Phys. 1977, 67, 4517-4533)
def get_ipc_score(mols):
	IPC_scores = []
	for i, mol in enumerate(mols):
		score = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				score += graph.Ipc(mol)
			except:
				pass

		IPC_scores.append(score)
	return IPC_scores

# Ertl SA_Score (J. Cheminform. 2009, 1, 8)
def get_sa_score(mols):
	# read fragment scores from file
	sa_score.readFragmentScores("fpscores")
	SA_Scores = []
	for i, mol in enumerate(mols):
		score = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				score += sa_score.calculateScore(mol)
			except:
				pass
		SA_Scores.append(score)
	return SA_Scores

# Boettcher Score (J. Chem. Inf. Model. 2016, 56, 3, 462–470)
def get_boettcher_score(smis):
	obConversion = openbabel.OBConversion()
	obConversion.SetInAndOutFormats("smi", "smi")
	bottch = boettcher.BottchScore("False")
	Boettcher_Scores = []
	for i, smi in enumerate(smis):
		score = 0
		for s in smi.split('.'):
			try:
				obmol = openbabel.OBMol()
				obConversion.ReadString(obmol, s)
				score += bottch.score(obmol)
			except:
				pass

		Boettcher_Scores.append(score)
	return Boettcher_Scores

#  SCScore (J. Chem. Inf. Model. 2018, 58, 2, 252)
def get_scscore(smis):
	model = standalone_model_numpy.SCScorer()
	resource_package = "molcomplex"
	resource_path = '/'.join(('models', 'full_reaxys_model_1024bool', 'model.ckpt-10654.as_numpy.json.gz'))
	json_path = pkg_resources.resource_filename(resource_package, resource_path)
	model.restore(json_path)
	#model.restore(os.path.join('.', 'models', 'full_reaxys_model_1024bool', 'model.ckpt-10654.as_numpy.json.gz'))
	SC_Scores = []
	for i, smi in enumerate(smis):
		scscore = 0
		for s in smi.split('.'):
			try:
				(s, score) = model.get_score_from_smi(s)
				scscore += score
			except:
				pass
		SC_Scores.append(scscore)
	return SC_Scores

#  SYBA Score (J. Cheminformatics 2020, 12, 35)
def get_sybascore(smis):
	syba = SybaClassifier()
	syba.fitDefaultScore()
	SYBA_Scores = []
	for i, smi in enumerate(smis):
		sybascore = 0
		for s in smi.split('.'):
			try:
				score = syba.predict(smi)
				sybascore += score
			except Exception as E:
				#print(E)
				pass
		SYBA_Scores.append(sybascore)
	return SYBA_Scores


def get_rucker_twc(mols):
	twc_scores = []
	for i, mol in enumerate(mols):
		score = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				score += rucker_twc.twc(mol)
			except:
				pass
		twc_scores.append(score)
	return twc_scores

def get_proudfoot_index(mols):
	pi_scores = []
	for i, mol in enumerate(mols):
		score = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				score += proudfoot.proudfoot_index(mol)[0]
			except:
				pass
		pi_scores.append(score)
	return pi_scores

def get_graph(mol):
	Chem.Kekulize(mol)
	atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
	am = Chem.GetAdjacencyMatrix(mol, useBO=True)
	for i,atom in enumerate(atoms):
		am[i,i] = atom
	G = nx.from_numpy_array(am)
	return G

def get_graph_edit_distance(mols):
	graph_edit_distance = []
	graph_edit_distance.append(0)
	graph_base = get_graph(mols[0])
	for i in range(1, len(mols)):
		score = 0
		try:
			graph = get_graph(mols[i])
			score = nx.graph_edit_distance(graph_base, graph)
			print(score)
		except:
			pass
		graph_edit_distance.append(score)
	return graph_edit_distance

''' Taken from Merck Paper descriptors '''

''' 1. Descriptor complex - AP and TT '''
def DESCRIPTORCOMPLEXITY_UNIQUEAP(mols):
	num_uniq_ap_list = []
	for i, mol in enumerate(mols):
		num_uniq_ap = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				mol_noHs = Chem.RemoveHs(mol)
				mol_ap_fp = Pairs.GetAtomPairFingerprint(mol_noHs)
				num_uniq_ap += len(mol_ap_fp.GetNonzeroElements())
			except:
				pass
		num_uniq_ap_list.append(num_uniq_ap)
	return num_uniq_ap_list

def DESCRIPTORCOMPLEXITY_UNIQUETT(mols):
	# need to check as the example from merck as 19 unique but actual unique is 13.
	num_uniq_tt_list = []
	for i, mol in enumerate(mols):
		num_uniq_tt = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				mol_noHs = Chem.RemoveHs(mol)
				mol_tt_fp = Torsions.GetTopologicalTorsionFingerprint(mol_noHs)
				num_uniq_tt += len(mol_tt_fp.GetNonzeroElements())
			except:
				pass
		num_uniq_tt_list.append(num_uniq_tt)
	return num_uniq_tt_list

def DESCRIPTORCOMPLEXITY_TOTALAP(mols):
	num_tot_AP_list = []
	for i, mol in enumerate(mols):
		num_tot_AP = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				mol_noHs = Chem.RemoveHs(mol)
				num_noHs = mol_noHs.GetNumAtoms()
				num_tot_AP += (num_noHs*(num_noHs - 1))/2
			except:
				pass
		num_tot_AP_list.append(num_tot_AP)
	return num_tot_AP_list

def DESCRIPTORCOMPLEXITY_TOTALTT(mols):
	num_tot_TT_list = []
	for i, mol in enumerate(mols):
		num_tot_TT = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				num_tot_TT += len(Chem.rdmolops.FindAllPathsOfLengthN(mol,3))
			except:
				pass
		num_tot_TT_list.append(num_tot_TT)
	return num_tot_TT_list


def DESCRIPTORCOMPLEXITY_APCOMPLEX(mols):
	apc_list = []
	num_uniq = DESCRIPTORCOMPLEXITY_UNIQUEAP(mols)
	num_tot = DESCRIPTORCOMPLEXITY_TOTALAP(mols)
	for i,j in zip(num_uniq,num_tot):
		if j == 0: apc = 0
		else: apc = i/j
		apc_list.append(apc)
	return apc_list

def DESCRIPTORCOMPLEXITY_TTCOMPLEX(mols):
	ttc_list = []
	num_uniq = DESCRIPTORCOMPLEXITY_UNIQUETT(mols)
	num_tot = DESCRIPTORCOMPLEXITY_TOTALTT(mols)
	for i,j in zip(num_uniq,num_tot):
		if j == 0: ttc = 0
		else: ttc = i/j
		ttc_list.append(ttc)
	return ttc_list

''' 2. MOE_2D '''


''' 3. SP3CARBONS '''

def SP3CARBONS_TOTALATOM_COUNT(mols):
	totatoms_list = []
	for i, mol in enumerate(mols):
		totatoms = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				totatoms += mol.GetNumAtoms()
			except:
				pass
		totatoms_list.append(totatoms)
	return totatoms_list

def SP3CARBONS_TOTALCARBON_COUNT(mols):
	totcar_list = []
	for i, mol in enumerate(mols):
		totcar =  0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				for atom in mol.GetAtoms():
					if atom.GetSymbol() == 'C':
						totcar += 1
			except:
				pass
		totcar_list.append(totcar)
	return totcar_list


def SP3CARBONS_CAR_ALLATOM_RATIO(mols):
	cratio_list = []
	for i, mol in enumerate(mols):
		cratio = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			numcar = SP3CARBONS_TOTALCARBON_COUNT([mol])
			numatoms  = SP3CARBONS_TOTALATOM_COUNT([mol])
			for i,j in zip(numcar,numatoms):
				if j == 0: cratio += 0
				else: cratio += i/j
		cratio_list.append(cratio)
	return cratio_list

######################

def SP3CARBONS_CHIRAL_ALLATOM_RATIO(mols):
	chiratio_list = []
	for i, mol in enumerate(mols):
		chiratio = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			numchiral = SP3CARBONS_CHIRAL_COUNT([mol])
			numatoms  = SP3CARBONS_TOTALATOM_COUNT([mol])
			for i,j in zip(numchiral,numatoms):
				if j == 0: chiratio += 0
				else: chiratio += i/j
		chiratio_list.append(chiratio)
	return chiratio_list

def SP3CARBONS_CHIRAL_ALLCARBON_RATIO(mols):
	chiratio_list = []
	for i, mol in enumerate(mols):
		chiratio = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			numchiral = SP3CARBONS_CHIRAL_COUNT([mol])
			numcarbons  = SP3CARBONS_TOTALCARBON_COUNT([mol])
			for i,j in zip(numchiral,numcarbons):
				if j == 0: chiratio += 0
				else: chiratio += i/j
		chiratio_list.append(chiratio)
	return chiratio_list



def SP3CARBONS_CHIRAL_COUNT(mols):
	c_chiral_list = []
	for i, mol in enumerate(mols):
		c_chiral = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				chiralcenters = Chem.FindMolChiralCenters(mol,includeUnassigned=True)
				for list in chiralcenters:
					if mol.GetAtomWithIdx(list[0]).GetSymbol() == 'C':
						c_chiral += 1
			except:
				pass
		c_chiral_list.append(c_chiral)
	return c_chiral_list

#########################

def SP3CARBONS_CSP2_COUNT(mols):
	csp2_list = []
	for i, mol in enumerate(mols):
		csp2 =  0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				for atom in mol.GetAtoms():
					if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP2:
						csp2 += 1
			except:
				pass
		csp2_list.append(csp2)
	return csp2_list

def SP3CARBONS_CSP2_ALLATOM_RATIO(mols):
	csp2ratio_list = []
	for i, mol in enumerate(mols):
		csp2ratio = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			numcsp2 = SP3CARBONS_CSP2_COUNT([mol])
			numatoms  = SP3CARBONS_TOTALATOM_COUNT([mol])
			for i,j in zip(numcsp2,numatoms):
				if j == 0: csp2ratio += 0
				else: csp2ratio += i/j
		csp2ratio_list.append(csp2ratio)
	return csp2ratio_list


def SP3CARBONS_CSP2_ALLCARBON_RATIO(mols):
	csp2ratio_list = []
	for i, mol in enumerate(mols):
		csp2ratio = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			numcsp2 = SP3CARBONS_CSP2_COUNT([mol])
			numcar  = SP3CARBONS_TOTALCARBON_COUNT([mol])
			for i,j in zip(numcsp2,numcar):
				if j == 0: csp2ratio += 0
				else: csp2ratio += i/j
		csp2ratio_list.append(csp2ratio)
	return csp2ratio_list

##########################

def SP3CARBONS_CSP3_COUNT(mols):
	csp3_list = []
	for i, mol in enumerate(mols):
		csp3 =  0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				for atom in mol.GetAtoms():
					if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP3:
						csp3 += 1
			except:
				pass
		csp3_list.append(csp3)
	return csp3_list

def SP3CARBONS_CSP3_ALLATOM_RATIO(mols):
	csp3ratio_list = []
	for i, mol in enumerate(mols):
		csp3ratio = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			numcsp3 = SP3CARBONS_CSP3_COUNT([mol])
			numatoms  = SP3CARBONS_TOTALATOM_COUNT([mol])
			for i,j in zip(numcsp3,numatoms):
				if j == 0: csp3ratio += 0
				else: csp3ratio += i/j
		csp3ratio_list.append(csp3ratio)
	return csp3ratio_list

def SP3CARBONS_CSP3_ALLCARBON_RATIO(mols):
	csp3ratio_list = []
	for i, mol in enumerate(mols):
		csp3ratio = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			numcsp3 = SP3CARBONS_CSP3_COUNT([mol])
			numcar  = SP3CARBONS_TOTALCARBON_COUNT([mol])
			for i,j in zip(numcsp3,numcar):
				if j == 0: csp3ratio += 0
				else: csp3ratio += i/j
		csp3ratio_list.append(csp3ratio)
	return csp3ratio_list

###########################

def SP3CARBONS_CSP_COUNT(mols):
	csp_list = []
	for i, mol in enumerate(mols):
		csp =  0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				for atom in mol.GetAtoms():
					if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP:
						csp += 1
			except:
				pass
		csp_list.append(csp)
	return csp_list

def SP3CARBONS_CSP_ALLATOM_RATIO(mols):
	cspratio_list = []
	for i, mol in enumerate(mols):
		cspratio = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			numcsp = SP3CARBONS_CSP_COUNT([mol])
			numatoms  = SP3CARBONS_TOTALATOM_COUNT([mol])
			for i,j in zip(numcsp,numatoms):
				if j == 0: cspratio += 0
				else: cspratio += i/j
		cspratio_list.append(cspratio)
	return cspratio_list

def SP3CARBONS_CSP_ALLCARBON_RATIO(mols):
	cspratio_list = []
	for i, mol in enumerate(mols):
		cspratio = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			numcsp = SP3CARBONS_CSP_COUNT([mol])
			numcar  = SP3CARBONS_TOTALCARBON_COUNT([mol])
			for i,j in zip(numcsp,numcar):
				if j == 0: cspratio += 0
				else: cspratio += i/j
		cspratio_list.append(cspratio)
	return cspratio_list

''' general rdkit descriptors  - rings information'''

def RINGINFO_NUM_ALI_CARBOCYCLE(mols):
	aliccycle_list = []
	for i, mol in enumerate(mols):
		aliccycle = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				aliccycle += Chem.rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)
			except:
				pass
		aliccycle_list.append(aliccycle)
	return aliccycle_list

def RINGINFO_NUM_ALI_HETEROCYCLE(mols):
	alihetcycle_list = []
	for i, mol in enumerate(mols):
		alihetcycle = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				alihetcycle += Chem.rdMolDescriptors.CalcNumAliphaticHeterocycles(mol)
			except:
				pass
		alihetcycle_list.append(alihetcycle)
	return alihetcycle_list

def RINGINFO_NUM_ALI_RINGS(mols):
	alirings_list = []
	for i, mol in enumerate(mols):
		alirings = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				alirings += Chem.rdMolDescriptors.CalcNumAliphaticRings(mol)
			except:
				pass
		alirings_list.append(alirings)
	return alirings_list


def RINGINFO_NUM_ARO_CARBOCYCLE(mols):
	aroccycle_list = []
	for i, mol in enumerate(mols):
		aroccycle = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				aroccycle += Chem.rdMolDescriptors.CalcNumAromaticCarbocycles(mol)
			except:
				pass
		aroccycle_list.append(aroccycle)
	return aroccycle_list

def RINGINFO_NUM_ARO_HETEROCYCLE(mols):
	arohetcycle_list = []
	for i, mol in enumerate(mols):
		arohetcycle = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				arohetcycle += Chem.rdMolDescriptors.CalcNumAromaticHeterocycles(mol)
			except:
				pass
		arohetcycle_list.append(arohetcycle)
	return arohetcycle_list

def RINGINFO_NUM_ARO_RINGS(mols):
	arorings_list = []
	for i, mol in enumerate(mols):
		arorings = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				arorings += Chem.rdMolDescriptors.CalcNumAromaticRings(mol)
			except:
				pass
		arorings_list.append(arorings)
	return arorings_list

def RINGINFO_NUM_BRIDGE_ATOMS(mols):
	bridge_list = []
	for i, mol in enumerate(mols):
		bridge = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				bridge += Chem.rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
			except:
				pass
		bridge_list.append(bridge)
	return bridge_list

def RINGINFO_NUM_SPIRO_ATOMS(mols):
	spiro_list = []
	for i, mol in enumerate(mols):
		smi = Chem.MolToSmiles(mol)
		spiro = 0
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				spiro += Chem.rdMolDescriptors.CalcNumSpiroAtoms(mol)
			except:
				pass
		spiro_list.append(spiro)
	return spiro_list

''' RDKit descriptors - present '''

def WIENER_INDEX(mols):
	wieneridx_list = []
	for i, mol in enumerate(mols):
		wieneridx = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				amat = Chem.GetDistanceMatrix(mol)
				num_atoms = mol.GetNumAtoms()
				for i in range(num_atoms):
					for j in range(i+1,num_atoms):
						wieneridx += amat[i][j]
			except:
				pass
		wieneridx_list.append(wieneridx)
	return wieneridx_list

def SMILES_3_2(smis):
	smi32_list = []
	for i, smi in enumerate(smis):
		smi32 = 0
		for s in smi.split('.'):
			try:
				smi32 += len(s)**(3/2)
			except:
				pass
		smi32_list.append(smi32)
	return smi32_list

''' Taken from pubchem descriptors '''

def PUBCHEM_XLOGP(mols):
	xlogp_list = []
	for i, mol in enumerate(mols):
		xlogp = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				xlogp += Chem.rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
			except:
				pass
		xlogp_list.append(xlogp)
	return xlogp_list

def PUBCHEM_TPSA(mols):
	tpsa_list = []
	for i, mol in enumerate(mols):
		tpsa = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				tpsa += Chem.rdMolDescriptors.CalcTPSA(mol)
			except:
				pass
		tpsa_list.append(tpsa)
	return tpsa_list

def PUBCHEM_H_BOND_DONOR_COUNT(mols):
	hbd_list = []
	for i, mol in enumerate(mols):
		hbd = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				hbd += Chem.rdMolDescriptors.CalcNumHBD(mol)
			except:
				pass
		hbd_list.append(hbd)
	return hbd_list

def PUBCHEM_H_BOND_ACCEPTOR_COUNT(mols):
	hba_list = []
	for i, mol in enumerate(mols):
		hba = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				hba += Chem.rdMolDescriptors.CalcNumHBA(mol)
			except:
				pass
		hba_list.append(hba)
	return hba_list

def PUBCHEM_ROTATABLE_BOND_COUNT(mols):
	rbond_list = []
	for i, mol in enumerate(mols):
		rbond = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				rbond += Chem.rdMolDescriptors.CalcNumRotatableBonds(mol,strict=1)
			except:
				pass
		rbond_list.append(rbond)
	return rbond_list

def PUBCHEM_HEAVY_ATOM_COUNT(mols):
	hac_list = []
	for i, mol in enumerate(mols):
		hac = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				hac += Chem.Lipinski.HeavyAtomCount(mol)
			except:
				pass
		hac_list.append(hac)
	return hac_list

def PUBCHEM_ATOM_STEREO_COUNT(mols):
	asc_list = []
	for i, mol in enumerate(mols):
		asc = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				asc += Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol)
			except:
				pass
		asc_list.append(asc)
	return asc_list

def PUBCHEM_DEFINED_ATOM_STEREO_COUNT(mols):
	dasc_list = []
	all_count = PUBCHEM_ATOM_STEREO_COUNT(mols)
	undefined_count = PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(mols)
	for i,j in zip(all_count,undefined_count):
		dasc_list.append(i-j)
	return dasc_list

def PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(mols):
	uasc_list = []
	for i, mol in enumerate(mols):
		uasc = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				uasc += Chem.rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)
			except:
				pass
		uasc_list.append(uasc)
	return uasc_list

def PUBCHEM_BOND_STEREO_COUNT(mols):
	return np.zeros(len(mols))

def PUBCHEM_DEFINED_BOND_STEREO_COUNT(mols):
	return np.zeros(len(mols))

def PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(mols):
	return np.zeros(len(mols))

def PUBCHEM_COVALENT_UNIT_COUNT(mols):
	return np.zeros(len(mols))

''' Mordred provided descriptors '''

def KAPPA_SHAPE_INDEX1(mols):
	"""Defined on page 9 of https://www.epa.gov/sites/production/files/2015-05/documents/moleculardescriptorsguide-v102.pdf"""
	ksi1 = KappaShapeIndex.KappaShapeIndex1()
	ksi_list = []
	for i, mol in enumerate(mols):
		ksi = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				ksi += ksi1(mol)
			except:
				pass
		ksi_list.append(ksi)
	return ksi_list

def KAPPA_SHAPE_INDEX2(mols):
	ksi2 = KappaShapeIndex.KappaShapeIndex2()
	ksi_list = []
	for i, mol in enumerate(mols):
		ksi = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				ksi += ksi2(mol)
			except:
				pass
		ksi_list.append(ksi)
	return ksi_list

def KAPPA_SHAPE_INDEX3(mols):
	ksi3 = KappaShapeIndex.KappaShapeIndex3()
	ksi_list = []
	for i, mol in enumerate(mols):
		ksi = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				ksi += ksi3(mol)
			except:
				pass
		ksi_list.append(ksi)
	return ksi_list

def MCGOWAN_VOLUME(mols):
	mv = McGowanVolume.McGowanVolume()
	mv_num_list = []
	for i, mol in enumerate(mols):
		mv_num = 0 
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				mv_num += mv(mol)
			except:
				pass
		mv_num_list.append(mv_num)
	return mv_num_list

def MOE_TYPE_Labute_ASA(mols):
	lasa = MoeType.LabuteASA()
	lasa_num_list = []
	for i, mol in enumerate(mols):
		lasa_num = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				lasa_num += lasa(mol)
			except:
				pass
		lasa_num_list.append(lasa_num)
	return lasa_num_list

def MOE_TYPE_PEOE_VSA(mols):
	pvsa = MoeType.PEOE_VSA()
	pvsa_val_list = []
	for i, mol in enumerate(mols):
		pvsa_val = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				pvsa_val += pvsa(mol)
			except:
				pass
		pvsa_val_list.append(pvsa_val)
	return pvsa_val_list

def MOE_TYPE_SMR_VSA(mols):
	svsa = MoeType.SMR_VSA()
	svsa_val_list = []
	for i, mol in enumerate(mols):
		svsa_val = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				svsa_val += svsa(mol)
			except:
				pass
		svsa_val_list.append(svsa_val)
	return svsa_val_list

def MOE_TYPE_SLOGP_VSA(mols):
	slpvsa = MoeType.SlogP_VSA()
	slpvsa_val_list = []
	for i, mol in enumerate(mols):
		slpvsa_val = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				slpvsa_val += slpvsa(mol)
			except:
				pass
		slpvsa_val_list.append(slpvsa_val)
	return slpvsa_val_list

def MOE_TYPE_ESTATE_VSA(mols):
	esvsa = MoeType.EState_VSA()
	esvsa_val_list = []
	for i, mol in enumerate(mols):
		esvsa_val = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				esvsa_val += esvsa(mol)
			except:
				pass
		esvsa_val_list.append(esvsa_val)
	return esvsa_val_list

def VDW_VOLUME_ABC(mols):
	vvabc = VdwVolumeABC.VdwVolumeABC()
	vvabc_val_list = []
	for i, mol in enumerate(mols):
		vvabc_val = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				vvabc_val += vvabc(mol)
			except:
				pass
		vvabc_val_list.append(vvabc_val)
	return vvabc_val_list

def ZAGREB_INDEX(mols):
	zi = ZagrebIndex.ZagrebIndex()
	zi_val_list = []
	for i, mol in enumerate(mols):
		zi_val = 0
		smi = Chem.MolToSmiles(mol)
		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
			try:
				zi_val += zi(mol)
			except:
				pass
		zi_val_list.append(zi_val)
	return zi_val_list

# ''' DBSTEP descriptors - need 3d Coords'''

# def MOL_VOLUME(mols):
# 	mol_vol_list = []
# 	for i, mol in enumerate(mols):
# 		smi = Chem.MolToSmiles(mol)
# 		for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
# 		try:
# 			sterics = db.dbstep(mol, commandline=True)
# 			mol_vol = sterics.occ_vol
# 		except:
# 			pass
# 		mol_vol_list.append(mol_vol)
# 	return mol_vol_list

#template for new descriptors
# _list = []
# for i, mol in enumerate(mols):
		# smi = Chem.MolToSmiles(mol)
		# for mol in [Chem.MolFromSmiles(s) for s in smi.split('.')]:
# 	try:
#
# 	except:
# 		pass
# 	_list.append()
# return _list
