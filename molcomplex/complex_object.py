# rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from pandas import DataFrame
from molcomplex.descriptors import *

import pandas as pd
pd.option_context('display.max_rows', None, 'display.max_columns', None)
import warnings
warnings.simplefilter(action='ignore', category=UserWarning) #ignore pandas UserWarnings for storing information in molcomplex objects

class mol_complex(DataFrame):
    #definition of mol_complex objects
    def __init__(self,mol_smiles, twc=False):
        if isinstance(mol_smiles, (DataFrame, pd.core.internals.BlockManager) ):
            super(mol_complex, self).__init__(mol_smiles)
        else:
            super(mol_complex,self).__init__(mol_smiles, columns=["SMILES"])

        self.twc=twc
        mol_smiles = list(mol_smiles)

        #general names for molecule
        self.mol_objects = [Chem.MolFromSmiles(smi) for smi in mol_smiles]

        self.canonical_smiles = [Chem.rdmolfiles.MolToSmiles(mol) for mol in self.mol_objects]
        # self.inchi = [Chem.inchi.MolToInchi(mol) for mol in self.mol_objects]
        # self.inchikey =  [Chem.inchi.MolToInchiKey(mol) for mol in self.mol_objects]

        #simple molecular properties
        self['molecular_formula'] = [rdMolDescriptors.CalcMolFormula(mol) for mol in self.mol_objects]
        self['molecular_weight'] = [Descriptors.ExactMolWt(mol) for mol in self.mol_objects]

        #implemented molecular scores
        self['BALABAN'] = get_balaban_score(self.mol_objects)
        self['BERTZ'] = get_bertz_score(self.mol_objects)
        self['BOETTCHER'] = get_boettcher_score(self.canonical_smiles)
        self['HKALPHA'] = get_hallkieralpha_score(self.mol_objects)
        self['IPC'] = get_ipc_score(self.mol_objects)
        self['SAS'] = get_sa_score(self.mol_objects)
        self['SCS'] = get_scscore(self.canonical_smiles)
        self['SYBA'] = get_sybascore(self.canonical_smiles)
        self['PI'] = get_proudfoot_index(self.mol_objects)
        self['SPATIALSCORE'] = get_spatial_score(self.canonical_smiles)
        if self.twc:
            self['R-TWC'] = get_rucker_twc(self.mol_objects)

        #assessing functions for mol Descriptors from descriptors.py
        self['DESCRIPTORCOMPLEXITY_UNIQUEAP'] = DESCRIPTORCOMPLEXITY_UNIQUEAP(self.mol_objects)
        self['DESCRIPTORCOMPLEXITY_UNIQUETT'] = DESCRIPTORCOMPLEXITY_UNIQUETT(self.mol_objects)
        self['DESCRIPTORCOMPLEXITY_TOTALAP'] = DESCRIPTORCOMPLEXITY_TOTALAP(self.mol_objects)
        self['DESCRIPTORCOMPLEXITY_TOTALTT'] = DESCRIPTORCOMPLEXITY_TOTALTT(self.mol_objects)
        self['DESCRIPTORCOMPLEXITY_APCOMPLEX'] = DESCRIPTORCOMPLEXITY_APCOMPLEX(self.mol_objects)
        self['DESCRIPTORCOMPLEXITY_TTCOMPLEX'] = DESCRIPTORCOMPLEXITY_TTCOMPLEX(self.mol_objects)
        self['SP3CARBONS_TOTALATOM_COUNT'] = SP3CARBONS_TOTALATOM_COUNT(self.mol_objects)
        self['SP3CARBONS_TOTALCARBON_COUNT'] = SP3CARBONS_TOTALCARBON_COUNT(self.mol_objects)
        self['SP3CARBONS_CAR_ALLATOM_RATIO'] = SP3CARBONS_CAR_ALLATOM_RATIO(self.mol_objects)
        self['SP3CARBONS_CHIRAL_ALLATOM_RATIO'] = SP3CARBONS_CHIRAL_ALLATOM_RATIO(self.mol_objects)
        self['SP3CARBONS_CHIRAL_ALLCARBON_RATIO'] = SP3CARBONS_CHIRAL_ALLCARBON_RATIO(self.mol_objects)
        self['SP3CARBONS_CHIRAL_COUNT'] = SP3CARBONS_CHIRAL_COUNT(self.mol_objects)
        self['SP3CARBONS_CSP2_COUNT'] = SP3CARBONS_CSP2_COUNT(self.mol_objects)
        self['SP3CARBONS_CSP2_ALLATOM_RATIO'] = SP3CARBONS_CSP2_ALLATOM_RATIO(self.mol_objects)
        self['SP3CARBONS_CSP2_ALLCARBON_RATIO'] = SP3CARBONS_CSP2_ALLCARBON_RATIO(self.mol_objects)
        self['SP3CARBONS_CSP3_COUNT'] = SP3CARBONS_CSP3_COUNT(self.mol_objects)
        self['SP3CARBONS_CSP3_ALLATOM_RATIO'] = SP3CARBONS_CSP3_ALLATOM_RATIO(self.mol_objects)
        self['SP3CARBONS_CSP3_ALLCARBON_RATIO'] = SP3CARBONS_CSP3_ALLCARBON_RATIO(self.mol_objects)
        self['SP3CARBONS_CSP_COUNT'] = SP3CARBONS_CSP_COUNT(self.mol_objects)
        self['SP3CARBONS_CSP_ALLATOM_RATIO'] = SP3CARBONS_CSP_ALLATOM_RATIO(self.mol_objects)
        self['SP3CARBONS_CSP_ALLCARBON_RATIO'] = SP3CARBONS_CSP_ALLCARBON_RATIO(self.mol_objects)
        self['RINGINFO_NUM_ALI_CARBOCYCLE'] = RINGINFO_NUM_ALI_CARBOCYCLE(self.mol_objects)
        self['RINGINFO_NUM_ALI_HETEROCYCLE'] = RINGINFO_NUM_ALI_HETEROCYCLE(self.mol_objects)
        self['RINGINFO_NUM_ALI_RINGS'] = RINGINFO_NUM_ALI_RINGS(self.mol_objects)
        self['RINGINFO_NUM_ARO_CARBOCYCLE'] = RINGINFO_NUM_ARO_CARBOCYCLE(self.mol_objects)
        self['RINGINFO_NUM_ARO_HETEROCYCLE'] = RINGINFO_NUM_ARO_HETEROCYCLE(self.mol_objects)
        self['RINGINFO_NUM_ARO_RINGS'] = RINGINFO_NUM_ARO_RINGS(self.mol_objects)
        self['RINGINFO_NUM_BRIDGE_ATOMS'] = RINGINFO_NUM_BRIDGE_ATOMS(self.mol_objects)
        self['RINGINFO_NUM_SPIRO_ATOMS'] = RINGINFO_NUM_SPIRO_ATOMS(self.mol_objects)
        self['WIENER_INDEX'] = WIENER_INDEX(self.mol_objects)
        self['SMILES_3_2'] = SMILES_3_2(self.canonical_smiles)
        self['PUBCHEM_XLOGP'] = PUBCHEM_XLOGP(self.mol_objects)
        self['PUBCHEM_TPSA'] = PUBCHEM_TPSA(self.mol_objects)
        self['PUBCHEM_H_BOND_DONOR_COUNT'] = PUBCHEM_H_BOND_DONOR_COUNT(self.mol_objects)
        self['PUBCHEM_H_BOND_ACCEPTOR_COUNT'] = PUBCHEM_H_BOND_ACCEPTOR_COUNT(self.mol_objects)
        self['PUBCHEM_ROTATABLE_BOND_COUNT'] = PUBCHEM_ROTATABLE_BOND_COUNT(self.mol_objects)
        self['PUBCHEM_HEAVY_ATOM_COUNT'] = PUBCHEM_HEAVY_ATOM_COUNT(self.mol_objects)
        self['PUBCHEM_ATOM_STEREO_COUNT'] = PUBCHEM_ATOM_STEREO_COUNT(self.mol_objects)
        self['PUBCHEM_DEFINED_ATOM_STEREO_COUNT'] = PUBCHEM_DEFINED_ATOM_STEREO_COUNT(self.mol_objects)
        self['PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT'] = PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(self.mol_objects)
        self['PUBCHEM_BOND_STEREO_COUNT'] = PUBCHEM_BOND_STEREO_COUNT(self.mol_objects)
        self['PUBCHEM_DEFINED_BOND_STEREO_COUNT'] = PUBCHEM_DEFINED_BOND_STEREO_COUNT(self.mol_objects)
        self['PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT'] = PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(self.mol_objects)
        self['PUBCHEM_COVALENT_UNIT_COUNT'] = PUBCHEM_COVALENT_UNIT_COUNT(self.mol_objects)
        self['KAPPA_SHAPE_INDEX1'] = KAPPA_SHAPE_INDEX1(self.mol_objects)
        self['KAPPA_SHAPE_INDEX2'] = KAPPA_SHAPE_INDEX2(self.mol_objects)
        self['KAPPA_SHAPE_INDEX3'] = KAPPA_SHAPE_INDEX3(self.mol_objects)
        self['MCGOWAN_VOLUME'] = MCGOWAN_VOLUME(self.mol_objects)
        self['MOE_TYPE_Labute_ASA'] = MOE_TYPE_Labute_ASA(self.mol_objects)
        self['MOE_TYPE_PEOE_VSA'] = MOE_TYPE_PEOE_VSA(self.mol_objects)
        self['MOE_TYPE_SMR_VSA'] = MOE_TYPE_SMR_VSA(self.mol_objects)
        self['MOE_TYPE_SLOGP_VSA'] = MOE_TYPE_SLOGP_VSA(self.mol_objects)
        self['MOE_TYPE_ESTATE_VSA'] = MOE_TYPE_ESTATE_VSA(self.mol_objects)
        self['VDW_VOLUME_ABC'] = VDW_VOLUME_ABC(self.mol_objects)
        self['ZAGREB_INDEX'] = ZAGREB_INDEX(self.mol_objects)