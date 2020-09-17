#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	        used for genrating a descriptors     	#
#####################################################.

import rdkit
from rdkit import Chem
from rdkit.Chem.AtomPairs import Pairs,Torsions
from openbabel import openbabel
from mordred import CPSA, KappaShapeIndex, McGowanVolume, MoeType, VdwVolumeABC, ZagrebIndex
import dbstep.Dbstep as db


''' Taken from Merck Paper descriptors '''

''' 1. Descriptor complex - AP and TT '''
def DESCRIPTORCOMPLEXITY_UNIQUEAP(mol):
    mol_noHs = Chem.RemoveHs(mol)
    mol_ap_fp = Pairs.GetAtomPairFingerprint(mol_noHs)
    num_uniq_ap = len(mol_ap_fp.GetNonzeroElements())
    return num_uniq_ap


def DESCRIPTORCOMPLEXITY_UNIQUETT(mol):
    # need ti check as the example from merck as 19 unique but actual unique is 13.
    mol_noHs = Chem.RemoveHs(mol)
    mol_tt_fp = Torsions.GetTopologicalTorsionFingerprint(mol_noHs)
    num_uniq_tt = len(mol_tt_fp.GetNonzeroElements())
    return num_uniq_tt

def DESCRIPTORCOMPLEXITY_TOTALAP(mol):
    mol_noHs = Chem.RemoveHs(mol)
    num_noHs = mol_noHs.GetNumAtoms()
    num_tot_AP = (num_noHs*(num_noHs - 1))/2
    return num_tot_AP

def DESCRIPTORCOMPLEXITY_TOTALTT(mol):
    #need to find a wways to get total TT
    return 0

def DESCRIPTORCOMPLEXITY_APCOMPLEX(mol):
    num_uniq = DESCRIPTORCOMPLEXITY_UNIQUEAP(mol)
    num_tot = DESCRIPTORCOMPLEXITY_TOTALAP(mol)
    if num_tot == 0: return 0
    else: return num_uniq/num_tot

def DESCRIPTORCOMPLEXITY_TTCOMPLEX(mol):
    num_uniq = DESCRIPTORCOMPLEXITY_UNIQUETT(mol)
    num_tot = DESCRIPTORCOMPLEXITY_TOTALTT(mol)
    if num_tot == 0: return 0
    else: return num_uniq/num_tot

''' 2. MOE_2D '''


''' 3. SP3CARBONS '''

def SP3CARBONS_TOTALATOM_COUNT(mol):
    totatoms  = mol.GetNumAtoms()
    return totatoms

def SP3CARBONS_TOTALCARBON_COUNT(mol):
    totcar =  0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            totcar += 1
    return totcar

def SP3CARBONS_CAR_ALLATOM_RATIO(mol):
    numcar = SP3CARBONS_TOTALCARBON_COUNT(mol)
    numatoms  = SP3CARBONS_TOTALATOM_COUNT(mol)
    return numcar/numatoms

######################

def SP3CARBONS_CHIRAL_ALLATOM_RATIO(mol):
    numchiral = SP3CARBONS_CHIRAL_COUNT(mol)
    numatoms  = SP3CARBONS_TOTALATOM_COUNT(mol)
    return numchiral/numatoms

def SP3CARBONS_CHIRAL_ALLCARBON_RATIO(mol):
    numchiral = SP3CARBONS_CHIRAL_COUNT(mol)
    numcarbons  = SP3CARBONS_TOTALCARBON_COUNT(mol)
    return numchiral/numcarbons

def SP3CARBONS_CHIRAL_COUNT(mol):
    c_chiral = 0
    chiralcenters = Chem.FindMolChiralCenters(mol,includeUnassigned=True)
    for list in chiralcenters:
        if mol.GetAtomWithIdx(list[0]).GetSymbol() == 'C':
            c_chiral += 1
    return c_chiral

#########################

def SP3CARBONS_CSP2_COUNT(mol):
    csp2 =  0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP2:
            car += 1
    return csp2

def SP3CARBONS_CSP2_ALLATOM_RATIO(mol):
    numcsp2 = SP3CARBONS_CSP2_COUNT(mol)
    numatoms  = SP3CARBONS_TOTALATOM_COUNT(mol)
    return numcsp2/numatoms

def SP3CARBONS_CSP2_ALLCARBON_RATIO(mol):
    numcsp2 = SP3CARBONS_CSP2_COUNT(mol)
    numcar  = SP3CARBONS_TOTALCARBON_COUNT(mol)
    return numcsp2/numcar

##########################

def SP3CARBONS_CSP3_COUNT(mol):
    csp3 =  0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP3:
            csp3 += 1
    return csp3

def SP3CARBONS_CSP3_ALLATOM_RATIO(mol):
    numcsp3 = SP3CARBONS_CSP3_COUNT(mol)
    numatoms  = SP3CARBONS_TOTALATOM_COUNT(mol)
    return numcsp3/numatoms

def SP3CARBONS_CSP3_ALLCARBON_RATIO(mol):
    numcsp3 = SP3CARBONS_CSP3_COUNT(mol)
    numcar  = SP3CARBONS_TOTALCARBON_COUNT(mol)
    return numcsp3/numcar

###########################

def SP3CARBONS_CSP_COUNT(mol):
    csp =  0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP:
            csp += 1
    return csp

def SP3CARBONS_CSP_ALLATOM_RATIO(mol):
    numcsp = SP3CARBONS_CSP_COUNT(mol)
    numatoms  = SP3CARBONS_TOTALATOM_COUNT(mol)
    return numcsp/numatoms

def SP3CARBONS_CSP_ALLCARBON_RATIO(mol):
    numcsp = SP3CARBONS_CSP_COUNT(mol)
    numcar  = SP3CARBONS_TOTALCARBON_COUNT(mol)
    return numcsp/numcar

''' general rdkit descriptors  - rings information'''

def RINGINFO_NUM_ALI_CARBOCYCLE(mol):
    return Chem.rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)

def RINGINFO_NUM_ALI_HETEROCYCLE(mol):
    return Chem.rdMolDescriptors.CalcNumAliphaticHeterocycles(mol)

def RINGINFO_NUM_ALI_RINGS(mol):
    return Chem.rdMolDescriptors.CalcNumAliphaticRings(mol)

def RINGINFO_NUM_ARO_CARBOCYCLE(mol):
    return Chem.rdMolDescriptors.CalcNumAromaticCarbocycles(mol)

def RINGINFO_NUM_ARO_HETEROCYCLE(mol):
    return Chem.rdMolDescriptors.CalcNumAromaticHeterocycles(mol)

def RINGINFO_NUM_ARO_RINGS(mol):
    return Chem.rdMolDescriptors.CalcNumAromaticRings(mol)

def RINGINFO_NUM_BRIDGE_ATOMS(mol):
    return Chem.rdMolDescriptors.CalcNumBridgeheadAtoms(mol)

def RINGINFO_NUM_SPIRO_ATOMS(mol):
    return Chem.rdMolDescriptors.CalcNumSpiroAtoms(mol)

''' RDKit descriptors - present '''

def WIENER_INDEX(mol):
    res = 0
    amat = Chem.GetDistanceMatrix(mol)
    num_atoms = mol.GetNumAtoms()
    for i in range(num_atoms):
        for j in range(i+1,num_atoms):
            res += amat[i][j]
    return res

def SMILES_3_2(smi):
    return len(smi)**(3/2)


''' Taken from pubchem descriptors '''

def PUBCHEM_XLOGP(mol):
    return Chem.rdMolDescriptors.CalcCrippenDescriptors(mol)[0]

def PUBCHEM_TPSA(mol):
    return Chem.rdMolDescriptors.CalcTPSA(mol)

def PUBCHEM_H_BOND_DONOR_COUNT(mol):
    return Chem.rdMolDescriptors.CalcNumHBD(mol)

def PUBCHEM_H_BOND_ACCEPTOR_COUNT(mol):
    return Chem.rdMolDescriptors.CalcNumHBA(mol)

def PUBCHEM_ROTATABLE_BOND_COUNT(mol):
    return Chem.rdMolDescriptors.CalcNumRotatableBonds(mol,strict=1)

def PUBCHEM_HEAVY_ATOM_COUNT(mol):
    return Chem.Lipinski.HeavyAtomCount(mol)

def PUBCHEM_ATOM_STEREO_COUNT(mol):
    return Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol)

def PUBCHEM_DEFINED_ATOM_STEREO_COUNT(mol):
    return PUBCHEM_ATOM_STEREO_COUNT(mol) - PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(mol)

def PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(mol):
    return Chem.rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)

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
    return ksi1(mol)
    
def KAPPA_SHAPE_INDEX2(mol):
    ksi2 = KappaShapeIndex.KappaShapeIndex2()
    return ksi2(mol)
    
def KAPPA_SHAPE_INDEX3(mol):
    ksi3 = KappaShapeIndex.KappaShapeIndex3()
    return ksi3(mol)
    
def MCGOWAN_VOLUME(mol):
    mv = McGowanVolume.McGowanVolume()
    return mv(mol)
    
def MOE_TYPE_Labute_ASA(mol):
    lasa = MoeType.LabuteASA()
    return lasa(mol)
    
def MOE_TYPE_PEOE_VSA(mol):
    pvsa = MoeType.PEOE_VSA()
    return pvsa(mol)
    
def MOE_TYPE_SMR_VSA(mol):
    svsa = MoeType.SMR_VSA()
    return svsa(mol)  

def MOE_TYPE_SLOGP_VSA(mol):
    slpvsa = MoeType.SlogP_VSA()
    return slpvsa(mol) 
    
def MOE_TYPE_ESTATE_VSA(mol):
    esvsa = MoeType.EState_VSA()
    return esvsa(mol) 

def VDW_VOLUME_ABC(mol):
    vvabc = VdwVolumeABC.VdwVolumeABC()
    return vvabc(mol)
    
def ZAGREB_INDEX(mol):
    zi = ZagrebIndex.ZagrebIndex()
    return zi(mol)

    
''' DBSTEP descriptors - need 3d Coords'''

def MOL_VOLUME(mol):
    sterics = db.dbstep(mol, commandline=True)
    return sterics.occ_vol
    
