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
#import dbstep.Dbstep as db

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

def get_rucker_twc(mols):
	twc_scores = []
	for i, mol in enumerate(mols):
		try:
			score = rucker_twc.twc(mol)
		except:
			 score = np.nan

		twc_scores.append(score)
	return twc_scores

def get_proudfoot_index(mols):
	pi_scores = []
	for i, mol in enumerate(mols):
		try:
			score = proudfoot.proudfoot_index(mol)[0]
		except:
			 score = np.nan
		pi_scores.append(score)
	return pi_scores

''' Taken from Merck Paper descriptors '''

''' 1. Descriptor complex - AP and TT '''
def DESCRIPTORCOMPLEXITY_UNIQUEAP(mols):
	num_uniq_ap_list = []
	for i, mol in enumerate(mols):
		try:
			mol_noHs = Chem.RemoveHs(mol)
			mol_ap_fp = Pairs.GetAtomPairFingerprint(mol_noHs)
			num_uniq_ap = len(mol_ap_fp.GetNonzeroElements())
		except:
			num_uniq_ap = np.nan
		num_uniq_ap_list.append(num_uniq_ap)
	return num_uniq_ap_list

def DESCRIPTORCOMPLEXITY_UNIQUETT(mols):
	# need to check as the example from merck as 19 unique but actual unique is 13.
	num_uniq_tt_list = []
	for i, mol in enumerate(mols):
		try:
			mol_noHs = Chem.RemoveHs(mol)
			mol_tt_fp = Torsions.GetTopologicalTorsionFingerprint(mol_noHs)
			num_uniq_tt = len(mol_tt_fp.GetNonzeroElements())
		except:
			num_uniq_tt = np.nan
		num_uniq_tt_list.append(num_uniq_tt)
	return num_uniq_tt_list

def DESCRIPTORCOMPLEXITY_TOTALAP(mols):
	num_tot_AP_list = []
	for i, mol in enumerate(mols):
		try:
			mol_noHs = Chem.RemoveHs(mol)
			num_noHs = mol_noHs.GetNumAtoms()
			num_tot_AP = (num_noHs*(num_noHs - 1))/2
		except:
			num_tot_AP= np.nan
		num_tot_AP_list.append(num_tot_AP)
	return num_tot_AP_list

def DESCRIPTORCOMPLEXITY_TOTALTT(mols):
	#need to find a wways to get total TT
	return 0

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
	for i, mol in enumerate(mols):
		try:
			num_uniq = DESCRIPTORCOMPLEXITY_UNIQUETT(mol)
			num_tot = DESCRIPTORCOMPLEXITY_TOTALTT(mol)
			if num_tot == 0: ttc = 0
			else: ttc = num_uniq/num_tot
		except:
			ttc = np.nan
		ttc_list.append(ttc)
	return ttc_list

''' 2. MOE_2D '''


''' 3. SP3CARBONS '''

def SP3CARBONS_TOTALATOM_COUNT(mols):
	totatoms_list = []
	for i, mol in enumerate(mols):
		try:
			totatoms = mol.GetNumAtoms()
		except:
			totatoms = np.nan
		totatoms_list.append(totatoms)
	return totatoms_list

def SP3CARBONS_TOTALCARBON_COUNT(mols):
	totcar_list = []
	for i, mol in enumerate(mols):
		try:
			totcar =  0
			for atom in mol.GetAtoms():
				if atom.GetSymbol() == 'C':
					totcar += 1
		except:
			totcar = np.nan
		totcar_list.append(totcar)
	return totcar_list

def SP3CARBONS_CAR_ALLATOM_RATIO(mols):
	cratio_list = []

	numcar = SP3CARBONS_TOTALCARBON_COUNT(mols)
	numatoms  = SP3CARBONS_TOTALATOM_COUNT(mols)
	for i,j in zip(numcar,numatoms):
		if j == 0: cratio = 0
		else: cratio = i/j
		cratio_list.append(cratio)
	return cratio_list

######################

def SP3CARBONS_CHIRAL_ALLATOM_RATIO(mols):
	chiratio_list = []
	numchiral = SP3CARBONS_CHIRAL_COUNT(mols)
	numatoms  = SP3CARBONS_TOTALATOM_COUNT(mols)
	for i,j in zip(numchiral,numatoms):
		if j == 0: chiratio = 0
		else: chiratio = i/j
		chiratio_list.append(chiratio)
	return chiratio_list

def SP3CARBONS_CHIRAL_ALLCARBON_RATIO(mols):
	chiratio_list = []
	numchiral = SP3CARBONS_CHIRAL_COUNT(mols)
	numcarbons  = SP3CARBONS_TOTALCARBON_COUNT(mols)
	for i,j in zip(numchiral,numcarbons):
		if j == 0: chiratio = 0
		else: chiratio = i/j
		chiratio_list.append(chiratio)
	return chiratio_list



def SP3CARBONS_CHIRAL_COUNT(mols):
	c_chiral_list = []
	for i, mol in enumerate(mols):
		try:
			c_chiral = 0
			chiralcenters = Chem.FindMolChiralCenters(mol,includeUnassigned=True)
			for list in chiralcenters:
				if mol.GetAtomWithIdx(list[0]).GetSymbol() == 'C':
					c_chiral += 1
		except:
			c_chiral = np.nan
		c_chiral_list.append(c_chiral)
	return c_chiral_list

#########################

def SP3CARBONS_CSP2_COUNT(mols):
	csp2_list = []
	for i, mol in enumerate(mols):
		try:
			csp2 =  0
			for atom in mol.GetAtoms():
				if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP2:
					csp2 += 1
		except:
			csp2 = np.nan
		csp2_list.append(csp2)
	return csp2_list

def SP3CARBONS_CSP2_ALLATOM_RATIO(mols):
	csp2ratio_list = []
	numcsp2 = SP3CARBONS_CSP2_COUNT(mols)
	numatoms  = SP3CARBONS_TOTALATOM_COUNT(mols)
	for i,j in zip(numcsp2,numatoms):
		if j == 0: csp2ratio = 0
		else: csp2ratio = i/j
		csp2ratio_list.append(csp2ratio)
	return csp2ratio_list


def SP3CARBONS_CSP2_ALLCARBON_RATIO(mols):
	csp2ratio_list = []
	numcsp2 = SP3CARBONS_CSP2_COUNT(mols)
	numcar  = SP3CARBONS_TOTALCARBON_COUNT(mols)
	for i,j in zip(numcsp2,numcar):
		if j == 0: csp2ratio = 0
		else: csp2ratio = i/j
		csp2ratio_list.append(csp2ratio)
	return csp2ratio_list

##########################

def SP3CARBONS_CSP3_COUNT(mols):
	csp3_list = []
	for i, mol in enumerate(mols):
		try:
			csp3 =  0
			for atom in mol.GetAtoms():
				if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP3:
					csp3 += 1
		except:
			csp3 = np.nan
		csp3_list.append(csp3)
	return csp3_list

def SP3CARBONS_CSP3_ALLATOM_RATIO(mols):
	csp3ratio_list = []
	numcsp3 = SP3CARBONS_CSP3_COUNT(mols)
	numatoms  = SP3CARBONS_TOTALATOM_COUNT(mols)
	for i,j in zip(numcsp3,numatoms):
		if j == 0: csp3ratio = 0
		else: csp3ratio = i/j
		csp3ratio_list.append(csp3ratio)
	return csp3ratio_list

def SP3CARBONS_CSP3_ALLCARBON_RATIO(mols):
	csp3ratio_list = []
	numcsp3 = SP3CARBONS_CSP3_COUNT(mols)
	numcar  = SP3CARBONS_TOTALCARBON_COUNT(mols)
	for i,j in zip(numcsp3,numcar):
		if j == 0: csp3ratio = 0
		else: csp3ratio = i/j
		csp3ratio_list.append(csp3ratio)
	return csp3ratio_list

###########################

def SP3CARBONS_CSP_COUNT(mols):
	csp_list = []
	for i, mol in enumerate(mols):
		try:
			csp =  0
			for atom in mol.GetAtoms():
				if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP:
					csp += 1
		except:
			csp = np.nan
		csp_list.append(csp)
	return csp_list

def SP3CARBONS_CSP_ALLATOM_RATIO(mols):
	cspratio_list = []
	numcsp = SP3CARBONS_CSP_COUNT(mols)
	numatoms  = SP3CARBONS_TOTALATOM_COUNT(mols)
	for i,j in zip(numcsp,numatoms):
		if j == 0: cspratio = 0
		else: cspratio = i/j
		cspratio_list.append(cspratio)
	return cspratio_list

def SP3CARBONS_CSP_ALLCARBON_RATIO(mols):
	cspratio_list = []
	numcsp = SP3CARBONS_CSP_COUNT(mols)
	numcar  = SP3CARBONS_TOTALCARBON_COUNT(mols)
	for i,j in zip(numcsp,numcar):
		if j == 0: cspratio = 0
		else: cspratio = i/j
		cspratio_list.append(cspratio)
	return cspratio_list

''' general rdkit descriptors  - rings information'''

def RINGINFO_NUM_ALI_CARBOCYCLE(mols):
	aliccycle_list = []
	for i, mol in enumerate(mols):
		try:
			aliccycle = Chem.rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)
		except:
			aliccycle = np.nan
		aliccycle_list.append(aliccycle)
	return aliccycle_list

def RINGINFO_NUM_ALI_HETEROCYCLE(mols):
	alihetcycle_list = []
	for i, mol in enumerate(mols):
		try:
			alihetcycle = Chem.rdMolDescriptors.CalcNumAliphaticHeterocycles(mol)
		except:
			alihetcycle = np.nan
		alihetcycle_list.append(alihetcycle)
	return alihetcycle_list

def RINGINFO_NUM_ALI_RINGS(mols):
	alirings_list = []
	for i, mol in enumerate(mols):
		try:
			alirings = Chem.rdMolDescriptors.CalcNumAliphaticRings(mol)
		except:
			alirings= np.nan
		alirings_list.append(alirings)
	return alirings_list


def RINGINFO_NUM_ARO_CARBOCYCLE(mols):
	aroccycle_list = []
	for i, mol in enumerate(mols):
		try:
			aroccycle = Chem.rdMolDescriptors.CalcNumAromaticCarbocycles(mol)
		except:
			aroccycle= np.nan
		aroccycle_list.append(aroccycle)
	return aroccycle_list

def RINGINFO_NUM_ARO_HETEROCYCLE(mols):
	arohetcycle_list = []
	for i, mol in enumerate(mols):
		try:
			arohetcycle = Chem.rdMolDescriptors.CalcNumAromaticHeterocycles(mol)
		except:
			arohetcycle= np.nan
		arohetcycle_list.append(arohetcycle)
	return arohetcycle_list

def RINGINFO_NUM_ARO_RINGS(mols):
	arorings_list = []
	for i, mol in enumerate(mols):
		try:
			arorings = Chem.rdMolDescriptors.CalcNumAromaticRings(mol)
		except:
			arorings= np.nan
		arorings_list.append(arorings)
	return arorings_list

def RINGINFO_NUM_BRIDGE_ATOMS(mols):
	bridge_list = []
	for i, mol in enumerate(mols):
		try:
			bridge = Chem.rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
		except:
			bridge = np.nan
		bridge_list.append(bridge)
	return bridge_list

def RINGINFO_NUM_SPIRO_ATOMS(mols):
	spiro_list = []
	for i, mol in enumerate(mols):
		try:
			spiro = Chem.rdMolDescriptors.CalcNumSpiroAtoms(mol)
		except:
			spiro= np.nan
		spiro_list.append(spiro)
	return spiro_list

''' RDKit descriptors - present '''

def WIENER_INDEX(mols):
	wieneridx_list = []
	for i, mol in enumerate(mols):
		try:
			wieneridx = 0
			amat = Chem.GetDistanceMatrix(mol)
			num_atoms = mol.GetNumAtoms()
			for i in range(num_atoms):
				for j in range(i+1,num_atoms):
					wieneridx += amat[i][j]
		except:
			wieneridx = np.nan
		wieneridx_list.append(wieneridx)
	return wieneridx_list

def SMILES_3_2(smis):
	smi32_list = []
	for i, smi in enumerate(smis):
		try:
			smi32 = len(smi)**(3/2)
		except:
			smi32 = np.nan
		smi32_list.append(smi32)
	return smi32_list

''' Taken from pubchem descriptors '''

def PUBCHEM_XLOGP(mols):
	xlogp_list = []
	for i, mol in enumerate(mols):
		try:
			xlogp = Chem.rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
		except:
			xlogp = np.nan
		xlogp_list.append(xlogp)
	return xlogp_list

def PUBCHEM_TPSA(mols):
	tpsa_list = []
	for i, mol in enumerate(mols):
		try:
			tpsa = Chem.rdMolDescriptors.CalcTPSA(mol)
		except:
			tpsa = np.nan
		tpsa_list.append(tpsa)
	return tpsa_list

def PUBCHEM_H_BOND_DONOR_COUNT(mols):
	hbd_list = []
	for i, mol in enumerate(mols):
		try:
			hbd = Chem.rdMolDescriptors.CalcNumHBD(mol)
		except:
			hbd = np.nan
		hbd_list.append(hbd)
	return hbd_list

def PUBCHEM_H_BOND_ACCEPTOR_COUNT(mols):
	hba_list = []
	for i, mol in enumerate(mols):
		try:
			hba = Chem.rdMolDescriptors.CalcNumHBA(mol)
		except:
			hba = np.nan
		hba_list.append(hba)
	return hba_list

def PUBCHEM_ROTATABLE_BOND_COUNT(mols):
	rbond_list = []
	for i, mol in enumerate(mols):
		try:
			rbond = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol,strict=1)
		except:
			rbond = np.nan
		rbond_list.append(rbond)
	return rbond_list

def PUBCHEM_HEAVY_ATOM_COUNT(mols):
	hac_list = []
	for i, mol in enumerate(mols):
		try:
			hac = Chem.Lipinski.HeavyAtomCount(mol)
		except:
			hac = np.nan
		hac_list.append(hac)
	return hac_list

def PUBCHEM_ATOM_STEREO_COUNT(mols):
	asc_list = []
	for i, mol in enumerate(mols):
		try:
			asc = Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol)
		except:
			asc = np.nan
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
		try:
			uasc = Chem.rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)
		except:
			uasc = np.nan
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
		try:
			ksi = ksi1(mol)
		except:
			ksi = np.nan
		ksi_list.append(ksi)
	return ksi_list

def KAPPA_SHAPE_INDEX2(mols):
	ksi2 = KappaShapeIndex.KappaShapeIndex2()
	ksi_list = []
	for i, mol in enumerate(mols):
		try:
			ksi = ksi2(mol)
		except:
			ksi = np.nan
		ksi_list.append(ksi)
	return ksi_list

def KAPPA_SHAPE_INDEX3(mols):
	ksi3 = KappaShapeIndex.KappaShapeIndex3()
	ksi_list = []
	for i, mol in enumerate(mols):
		try:
			ksi = ksi3(mol)
		except:
			ksi = np.nan
		ksi_list.append(ksi)
	return ksi_list

def MCGOWAN_VOLUME(mols):
	mv = McGowanVolume.McGowanVolume()
	mv_num_list = []
	for i, mol in enumerate(mols):
		try:
			mv_num = mv(mol)
		except:
			mv_num = np.nan
		mv_num_list.append(mv_num)
	return mv_num_list

def MOE_TYPE_Labute_ASA(mols):
	lasa = MoeType.LabuteASA()
	lasa_num_list = []
	for i, mol in enumerate(mols):
		try:
			lasa_num = lasa(mol)
		except:
			lasa_num= np.nan
		lasa_num_list.append(lasa_num)
	return lasa_num_list

def MOE_TYPE_PEOE_VSA(mols):
	pvsa = MoeType.PEOE_VSA()
	pvsa_val_list = []
	for i, mol in enumerate(mols):
		try:
			pvsa_val = pvsa(mol)
		except:
			pvsa_val = np.nan
		pvsa_val_list.append(pvsa_val)
	return pvsa_val_list

def MOE_TYPE_SMR_VSA(mols):
	svsa = MoeType.SMR_VSA()
	svsa_val_list = []
	for i, mol in enumerate(mols):
		try:
			svsa_val =  svsa(mol)
		except:
			svsa_val= np.nan
		svsa_val_list.append(svsa_val)
	return svsa_val_list

def MOE_TYPE_SLOGP_VSA(mols):
	slpvsa = MoeType.SlogP_VSA()
	slpvsa_val_list = []
	for i, mol in enumerate(mols):
		try:
			slpvsa_val = slpvsa(mol)
		except:
			slpvsa_val= np.nan
		slpvsa_val_list.append(slpvsa_val)
	return slpvsa_val_list

def MOE_TYPE_ESTATE_VSA(mols):
	esvsa = MoeType.EState_VSA()
	esvsa_val_list = []
	for i, mol in enumerate(mols):
		try:
			esvsa_val = esvsa(mol)
		except:
			esvsa_val= np.nan
		esvsa_val_list.append(esvsa_val)
	return esvsa_val_list

def VDW_VOLUME_ABC(mols):
	vvabc = VdwVolumeABC.VdwVolumeABC()
	vvabc_val_list = []
	for i, mol in enumerate(mols):
		try:
			vvabc_val = vvabc(mol)
		except:
			vvabc_val= np.nan
		vvabc_val_list.append(vvabc_val)
	return vvabc_val_list

def ZAGREB_INDEX(mols):
	zi = ZagrebIndex.ZagrebIndex()
	zi_val_list = []
	for i, mol in enumerate(mols):
		try:
			zi_val = zi(mol)
		except:
			zi_val= np.nan
		zi_val_list.append(zi_val)
	return zi_val_list

''' DBSTEP descriptors - need 3d Coords'''

def MOL_VOLUME(mols):
	mol_vol_list = []
	for i, mol in enumerate(mols):
		try:
			sterics = db.dbstep(mol, commandline=True)
			mol_vol = sterics.occ_vol
		except:
			mol_vol	= np.nan
		mol_vol_list.append(mol_vol)
	return mol_vol_list

#template for new descriptors
# _list = []
# for i, mol in enumerate(mols):
# 	try:
#
# 	except:
# 		= np.nan
# 	_list.append()
# return _list
