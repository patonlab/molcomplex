import rdkit
from rdkit import Chem
from openbabel import openbabel


''' Taken from Merck Paper descriptors '''

''' 1. Descriptor complex - AP and TT '''
def DESCRIPTORCOMPLEXITY_APCOMPLEX(mol):
    return 0

def DESCRIPTORCOMPLEXITY_TTCOMPLEX(mol):
    return 0

def DESCRIPTORCOMPLEXITY_UNIQUEAP(mol):
    return 0

def DESCRIPTORCOMPLEXITY_UNIQUETT(mol):
    return 0

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
    numcarbons  = SP3CARBONS_CAR_COUNT(mol)
    return numchiral/numcarbons

def SP3CARBONS_CHIRAL_COUNT(mol):
    c_chiral = 0
    chiralcenters = Chem.FindMolChiralCenters(mol)
    for list in chiralcenters:
        if mol.GetAtomWithIdx(list[0]).GetSymbol() = 'C':
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


''' Taken from pubchem descriptors '''

def XLOGP(mol):
    return Chem.rdMolDescriptors.CalcCrippenDescriptors(mol)[0]

def TPSA(mol):
    return Chem.rdMolDescriptors.CalcTPSA(mol)

def H_BOND_DONOR_COUNT(mol):
    return Chem.rdMolDescriptors.CalcNumHBD(mol)

def H_BOND_ACCEPTOR_COUNT(mol):
    return Chem.rdMolDescriptors.CalcNumHBA(mol)

def ROTATABLE_BOND_COUNT(mol):
    return Chem.rdMolDescriptors.CalcNumRotatableBonds(mol,strict=1)

def HEAVY_ATOM_COUNT(mol):
    return Chem.Lipinski.HeavyAtomCount(mol)

def ATOM_STEREO_COUNT(mol):
    return Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol)

def DEFINED_ATOM_STEREO_COUNT(mol):
    return ATOM_STEREO_COUNT(mol) - UNDEFINED_ATOM_STEREO_COUNT(mol)

def UNDEFINED_ATOM_STEREO_COUNT(mol):
    return Chem.rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)

def BOND_STEREO_COUNT(mol):
    return 0

def DEFINED_BOND_STEREO_COUNT(mol):
    return 0

def UNDEFINED_ATOM_STEREO_COUNT(mol):
    return 0

def COVALENT_UNIT_COUNT(mol):
    return 0
