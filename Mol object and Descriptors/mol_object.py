# rdkit
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from pandas import DataFrame
from descriptors import *
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

#import the descriptors functions

class mol_complex(DataFrame):
    #definition of molobjects
    def __init__(self,smi):
        super().__init__([{"SMILES":smi}])
        #self.data = {}
        #general names for molecule
        #self.input_smiles = smi
        self.mol = Chem.MolFromSmiles(smi)
        
        self.canonical_smiles = Chem.rdmolfiles.MolToSmiles(self.mol)
        self.inchi = Chem.inchi.MolToInchi(self.mol)
        self.inchikey =  Chem.inchi.MolToInchiKey(self.mol)
        # self.iupac_name
        
        #simple molecular properties
        self['molecular_formula'] = rdMolDescriptors.CalcMolFormula(self.mol)
        self['molecular_weight'] = Descriptors.ExactMolWt(self.mol)
        
        #assessing functions for mol Descriptors from descriptors.py
        self['DESCRIPTORCOMPLEXITY_UNIQUEAP'] = DESCRIPTORCOMPLEXITY_UNIQUEAP(self.mol)
        self['DESCRIPTORCOMPLEXITY_UNIQUETT'] = DESCRIPTORCOMPLEXITY_UNIQUETT(self.mol)
        self['DESCRIPTORCOMPLEXITY_TOTALAP'] = DESCRIPTORCOMPLEXITY_TOTALAP(self.mol)
        self['DESCRIPTORCOMPLEXITY_TOTALTT'] = DESCRIPTORCOMPLEXITY_TOTALTT(self.mol)
        self['DESCRIPTORCOMPLEXITY_APCOMPLEX'] = DESCRIPTORCOMPLEXITY_APCOMPLEX(self.mol)
        self['DESCRIPTORCOMPLEXITY_TTCOMPLEX'] = DESCRIPTORCOMPLEXITY_TTCOMPLEX(self.mol)
        self['SP3CARBONS_TOTALATOM_COUNT'] = SP3CARBONS_TOTALATOM_COUNT(self.mol)
        self['SP3CARBONS_TOTALCARBON_COUNT'] = SP3CARBONS_TOTALCARBON_COUNT(self.mol)
        self['SP3CARBONS_CAR_ALLATOM_RATIO'] = SP3CARBONS_CAR_ALLATOM_RATIO(self.mol)
        self['SP3CARBONS_CHIRAL_ALLATOM_RATIO'] = SP3CARBONS_CHIRAL_ALLATOM_RATIO(self.mol)
        self['SP3CARBONS_CHIRAL_ALLCARBON_RATIO'] = SP3CARBONS_CHIRAL_ALLCARBON_RATIO(self.mol)
        self['SP3CARBONS_CHIRAL_COUNT'] = SP3CARBONS_CHIRAL_COUNT(self.mol)
        self['SP3CARBONS_CSP2_COUNT'] = SP3CARBONS_CSP2_COUNT(self.mol)
        self['SP3CARBONS_CSP2_ALLATOM_RATIO'] = SP3CARBONS_CSP2_ALLATOM_RATIO(self.mol)
        self['SP3CARBONS_CSP2_ALLCARBON_RATIO'] = SP3CARBONS_CSP2_ALLCARBON_RATIO(self.mol)
        self['SP3CARBONS_CSP3_COUNT'] = SP3CARBONS_CSP3_COUNT(self.mol)
        self['SP3CARBONS_CSP3_ALLATOM_RATIO'] = SP3CARBONS_CSP3_ALLATOM_RATIO(self.mol)
        self['SP3CARBONS_CSP3_ALLCARBON_RATIO'] = SP3CARBONS_CSP3_ALLCARBON_RATIO(self.mol)
        self['SP3CARBONS_CSP_COUNT'] = SP3CARBONS_CSP_COUNT(self.mol)
        self['SP3CARBONS_CSP_ALLATOM_RATIO'] = SP3CARBONS_CSP_ALLATOM_RATIO(self.mol)
        self['SP3CARBONS_CSP_ALLCARBON_RATIO'] = SP3CARBONS_CSP_ALLCARBON_RATIO(self.mol)
        self['RINGINFO_NUM_ALI_CARBOCYCLE'] = RINGINFO_NUM_ALI_CARBOCYCLE(self.mol)
        self['RINGINFO_NUM_ALI_HETEROCYCLE'] = RINGINFO_NUM_ALI_HETEROCYCLE(self.mol)
        self['RINGINFO_NUM_ALI_RINGS'] = RINGINFO_NUM_ALI_RINGS(self.mol)
        self['RINGINFO_NUM_ARO_CARBOCYCLE'] = RINGINFO_NUM_ARO_CARBOCYCLE(self.mol)
        self['RINGINFO_NUM_ARO_HETEROCYCLE'] = RINGINFO_NUM_ARO_HETEROCYCLE(self.mol)
        self['RINGINFO_NUM_ARO_RINGS'] = RINGINFO_NUM_ARO_RINGS(self.mol)
        self['RINGINFO_NUM_BRIDGE_ATOMS'] = RINGINFO_NUM_BRIDGE_ATOMS(self.mol)
        self['RINGINFO_NUM_SPIRO_ATOMS'] = RINGINFO_NUM_SPIRO_ATOMS(self.mol)
        self['WIENER_INDEX'] = WIENER_INDEX(self.mol)
        self['SMILES_3_2'] = SMILES_3_2(smi)
        self['PUBCHEM_XLOGP'] = PUBCHEM_XLOGP(self.mol)
        self['PUBCHEM_TPSA'] = PUBCHEM_TPSA(self.mol)
        self['PUBCHEM_H_BOND_DONOR_COUNT'] = PUBCHEM_H_BOND_DONOR_COUNT(self.mol)
        self['PUBCHEM_H_BOND_ACCEPTOR_COUNT'] = PUBCHEM_H_BOND_ACCEPTOR_COUNT(self.mol)
        self['PUBCHEM_ROTATABLE_BOND_COUNT'] = PUBCHEM_ROTATABLE_BOND_COUNT(self.mol)
        self['PUBCHEM_HEAVY_ATOM_COUNT'] = PUBCHEM_HEAVY_ATOM_COUNT(self.mol)
        self['PUBCHEM_ATOM_STEREO_COUNT'] = PUBCHEM_ATOM_STEREO_COUNT(self.mol)
        self['PUBCHEM_DEFINED_ATOM_STEREO_COUNT'] = PUBCHEM_DEFINED_ATOM_STEREO_COUNT(self.mol)
        self['PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT'] = PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(self.mol)
        self['PUBCHEM_BOND_STEREO_COUNT'] = PUBCHEM_BOND_STEREO_COUNT(self.mol)
        self['PUBCHEM_DEFINED_BOND_STEREO_COUNT'] = PUBCHEM_DEFINED_BOND_STEREO_COUNT(self.mol)
        self['PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT'] = PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(self.mol)
        self['PUBCHEM_COVALENT_UNIT_COUNT'] = PUBCHEM_COVALENT_UNIT_COUNT(self.mol)
        self['KAPPA_SHAPE_INDEX1'] = KAPPA_SHAPE_INDEX1(self.mol)
        self['KAPPA_SHAPE_INDEX2'] = KAPPA_SHAPE_INDEX2(self.mol)
        self['KAPPA_SHAPE_INDEX3'] = KAPPA_SHAPE_INDEX3(self.mol)
        self['MCGOWAN_VOLUME'] = MCGOWAN_VOLUME(self.mol)
        self['MOE_TYPE_Labute_ASA'] = MOE_TYPE_Labute_ASA(self.mol)
        self['MOE_TYPE_PEOE_VSA'] = MOE_TYPE_PEOE_VSA(self.mol)
        self['MOE_TYPE_SMR_VSA'] = MOE_TYPE_SMR_VSA(self.mol)
        self['MOE_TYPE_SLOGP_VSA'] = MOE_TYPE_SLOGP_VSA(self.mol)
        self['MOE_TYPE_ESTATE_VSA'] = MOE_TYPE_ESTATE_VSA(self.mol)
        self['VDW_VOLUME_ABC'] = VDW_VOLUME_ABC(self.mol)
        self['ZAGREB_INDEX'] = ZAGREB_INDEX(self.mol)

        
    # this method is makes it so our methods return an instance
    # of a mol object, instead of a regular DataFrame
    # obtained from https://dev.to/pj_trainor/extending-the-pandas-dataframe-133l
    @property
    def _constructor(self):
        return mol_complex

if __name__ == "__main__":

    #just for testing the class is working
    smi = 'CCCC'
    dict={'SMILES':smi}
    mol_object = mol_complex(smi)
    print('SMILES',mol_object.SMILES[0])
    print('MOLFromSMILES',mol_object.mol)
    print(mol_object)
    mol_object.to_csv("test_butane.csv",index=False)
