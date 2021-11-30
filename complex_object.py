# rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from pandas import DataFrame
from descriptors import *
try:
	from .metrics import standalone_model_numpy
except:
	from metrics import standalone_model_numpy
import pandas as pd
pd.option_context('display.max_rows', None, 'display.max_columns', None)
import warnings
warnings.simplefilter(action='ignore', category=UserWarning) #ignore pandas UserWarnings for storing information in molcomplex objects


class mol_complex(DataFrame):
    func_dict = {'BALABAN': get_balaban_score, 'BERTZ': get_bertz_score, 
        'BOETTCHER': get_boettcher_score, 'HKALPHA': get_hallkieralpha_score,
        'IPC': get_ipc_score, 'SAS': get_sa_score, 'SCS': get_scscore,
        'PI': get_proudfoot_index, 
        'DESCRIPTORCOMPLEXITY_UNIQUEAP': DESCRIPTORCOMPLEXITY_UNIQUEAP,
        'DESCRIPTORCOMPLEXITY_UNIQUETT': DESCRIPTORCOMPLEXITY_UNIQUETT,
        'DESCRIPTORCOMPLEXITY_TOTALAP': DESCRIPTORCOMPLEXITY_TOTALAP,
        'DESCRIPTORCOMPLEXITY_TOTALTT': DESCRIPTORCOMPLEXITY_TOTALTT,
        'DESCRIPTORCOMPLEXITY_APCOMPLEX': DESCRIPTORCOMPLEXITY_APCOMPLEX,
        'DESCRIPTORCOMPLEXITY_TTCOMPLEX': DESCRIPTORCOMPLEXITY_TTCOMPLEX,
        'SP3CARBONS_TOTALATOM_COUNT': SP3CARBONS_TOTALATOM_COUNT,
        'SP3CARBONS_TOTALCARBON_COUNT': SP3CARBONS_TOTALCARBON_COUNT,
        'SP3CARBONS_CAR_ALLATOM_RATIO': SP3CARBONS_CAR_ALLATOM_RATIO,
        'SP3CARBONS_CHIRAL_ALLATOM_RATIO': SP3CARBONS_CHIRAL_ALLATOM_RATIO,
        'SP3CARBONS_CHIRAL_ALLCARBON_RATIO': SP3CARBONS_CHIRAL_ALLCARBON_RATIO,
        'SP3CARBONS_CHIRAL_COUNT': SP3CARBONS_CHIRAL_COUNT,
        'SP3CARBONS_CSP2_COUNT': SP3CARBONS_CSP2_COUNT,
        'SP3CARBONS_CSP2_ALLATOM_RATIO': SP3CARBONS_CSP2_ALLATOM_RATIO,
        'SP3CARBONS_CSP2_ALLCARBON_RATIO': SP3CARBONS_CSP2_ALLCARBON_RATIO,
        'SP3CARBONS_CSP3_COUNT': SP3CARBONS_CSP3_COUNT,
        'SP3CARBONS_CSP3_ALLATOM_RATIO': SP3CARBONS_CSP3_ALLATOM_RATIO,
        'SP3CARBONS_CSP3_ALLCARBON_RATIO': SP3CARBONS_CSP3_ALLCARBON_RATIO,
        'SP3CARBONS_CSP_COUNT': SP3CARBONS_CSP_COUNT,
        'SP3CARBONS_CSP_ALLATOM_RATIO': SP3CARBONS_CSP_ALLATOM_RATIO,
        'SP3CARBONS_CSP_ALLCARBON_RATIO': SP3CARBONS_CSP_ALLCARBON_RATIO,
        'RINGINFO_NUM_ALI_CARBOCYCLE': RINGINFO_NUM_ALI_CARBOCYCLE,
        'RINGINFO_NUM_ALI_HETEROCYCLE': RINGINFO_NUM_ALI_HETEROCYCLE,
        'RINGINFO_NUM_ALI_RINGS': RINGINFO_NUM_ALI_RINGS,
        'RINGINFO_NUM_ARO_CARBOCYCLE': RINGINFO_NUM_ARO_CARBOCYCLE,
        'RINGINFO_NUM_ARO_HETEROCYCLE': RINGINFO_NUM_ARO_HETEROCYCLE,
        'RINGINFO_NUM_ARO_RINGS': RINGINFO_NUM_ARO_RINGS,
        'RINGINFO_NUM_BRIDGE_ATOMS': RINGINFO_NUM_BRIDGE_ATOMS,
        'RINGINFO_NUM_SPIRO_ATOMS': RINGINFO_NUM_SPIRO_ATOMS,
        'WIENER_INDEX': WIENER_INDEX, 'SMILES_3_2': SMILES_3_2,
        'PUBCHEM_XLOGP': PUBCHEM_XLOGP, 'PUBCHEM_TPSA': PUBCHEM_TPSA,
        'PUBCHEM_H_BOND_DONOR_COUNT': PUBCHEM_H_BOND_DONOR_COUNT,
        'PUBCHEM_H_BOND_ACCEPTOR_COUNT': PUBCHEM_H_BOND_ACCEPTOR_COUNT,
        'PUBCHEM_ROTATABLE_BOND_COUNT': PUBCHEM_ROTATABLE_BOND_COUNT,
        'PUBCHEM_HEAVY_ATOM_COUNT': PUBCHEM_HEAVY_ATOM_COUNT,
        'PUBCHEM_ATOM_STEREO_COUNT': PUBCHEM_ATOM_STEREO_COUNT,
        'PUBCHEM_DEFINED_ATOM_STEREO_COUNT': PUBCHEM_DEFINED_ATOM_STEREO_COUNT,
        'PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT': PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT,
        'PUBCHEM_BOND_STEREO_COUNT': PUBCHEM_BOND_STEREO_COUNT,
        'PUBCHEM_DEFINED_BOND_STEREO_COUNT': PUBCHEM_DEFINED_BOND_STEREO_COUNT,
        'PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT': PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT,
        'PUBCHEM_COVALENT_UNIT_COUNT': PUBCHEM_COVALENT_UNIT_COUNT,
        'KAPPA_SHAPE_INDEX1': KAPPA_SHAPE_INDEX1,
        'KAPPA_SHAPE_INDEX2': KAPPA_SHAPE_INDEX2,
        'KAPPA_SHAPE_INDEX3': KAPPA_SHAPE_INDEX3,
        'MCGOWAN_VOLUME': MCGOWAN_VOLUME,
        'MOE_TYPE_Labute_ASA': MOE_TYPE_Labute_ASA,
        'MOE_TYPE_PEOE_VSA': MOE_TYPE_PEOE_VSA,
        'MOE_TYPE_SMR_VSA': MOE_TYPE_SMR_VSA,
        'MOE_TYPE_SLOGP_VSA': MOE_TYPE_SLOGP_VSA,
        'MOE_TYPE_ESTATE_VSA': MOE_TYPE_ESTATE_VSA,
        'VDW_VOLUME_ABC': VDW_VOLUME_ABC,
        'ZAGREB_INDEX': ZAGREB_INDEX,
        }  

    #definition of mol_complex objects
    def __init__(self,mol_smiles,options):
        if isinstance(mol_smiles, (DataFrame, pd.core.internals.BlockManager) ):
            super(mol_complex, self).__init__(mol_smiles)
        else:
            super(mol_complex,self).__init__(mol_smiles, columns=["SMILES"])
        self.options = options
        self.mol_smiles = list(mol_smiles)
        #self.inchi = [Chem.inchi.MolToInchi(mol) for mol in self.mol_objects]
        #self.inchikey =  [Chem.inchi.MolToInchiKey(mol) for mol in self.mol_objects]
        #self.iupac_name

        #general names for molecule
        self.populate_df()

    
    def populate_df(self):
        mol_smiles = self.mol_smiles
        mol_objects = [Chem.MolFromSmiles(smi) for smi in self.mol_smiles]   
        
        func_dict = self.func_dict
        #simple molecular properties
        self['molecular_formula'] = [rdMolDescriptors.CalcMolFormula(mol) for mol in mol_objects]
        self['molecular_weight'] = [Descriptors.ExactMolWt(mol) for mol in mol_objects]
        #canonical_smiles = [Chem.rdmolfiles.MolToSmiles(mol) for mol in mol_objects]
        
        # self = self.reindex(columns = self.columns.tolist()
        #                           + list(func_dict.keys()))

        for key in func_dict.keys():
            self[key] = 0
        if self.options.twc:
            self['R-TWC'] = 0
        print(self)

        # only load model once for speed improvement
        model = standalone_model_numpy.SCScorer()
        model.restore(os.path.join('.', 'models', 'full_reaxys_model_1024bool', 'model.ckpt-10654.as_numpy.json.gz'))

        for i,smi in enumerate(mol_smiles):
            #get score for each mol, sum contribs if '.' present
            for key in func_dict.keys():
                for s in smi.split('.'):
                    mol = Chem.MolFromSmiles(s)
                    
                    if key == 'SCS':
                        #print(key,func_dict[key](mol,model),type(func_dict[key](mol,model)))
                        self.loc[i,key] += func_dict[key](mol,model)
                    else:
                        #print(key,func_dict[key](mol),type(func_dict[key](mol)))
                        self.loc[i,key] += func_dict[key](mol)
            if self.options.twc:
                for s in smi.split('.'):
                    mol = Chem.MolFromSmiles(s)
                    self.loc[i,'R-TWC'] += get_rucker_twc(mol)

        print(self)


        #implemented molecular scores
        #self['BALABAN'] = get_balaban_score(self.mol_objects)
        #self['BERTZ'] = get_bertz_score(self.mol_objects)
        #self['BOETTCHER'] = get_boettcher_score(mol_smiles)
        #self['HKALPHA'] = get_hallkieralpha_score(self.mol_objects)
        #self['IPC'] = get_ipc_score(self.mol_objects)
        #self['SAS'] = get_sa_score(self.mol_objects)
        #self['SCS'] = get_scscore(mol_smiles)
        #self['PI'] = get_proudfoot_index(self.mol_objects)
        # if options.twc:
        #     self['R-TWC'] = get_rucker_twc(self.mol_objects)



        #assessing functions for mol Descriptors from descriptors.py
        #self['DESCRIPTORCOMPLEXITY_UNIQUEAP'] = DESCRIPTORCOMPLEXITY_UNIQUEAP(self.mol_objects)
        #self['DESCRIPTORCOMPLEXITY_UNIQUETT'] = DESCRIPTORCOMPLEXITY_UNIQUETT(self.mol_objects)
        #self['DESCRIPTORCOMPLEXITY_TOTALAP'] = DESCRIPTORCOMPLEXITY_TOTALAP(self.mol_objects)
        #self['DESCRIPTORCOMPLEXITY_TOTALTT'] = DESCRIPTORCOMPLEXITY_TOTALTT(self.mol_objects)
        #self['DESCRIPTORCOMPLEXITY_APCOMPLEX'] = DESCRIPTORCOMPLEXITY_APCOMPLEX(self.mol_objects)
        #self['DESCRIPTORCOMPLEXITY_TTCOMPLEX'] = DESCRIPTORCOMPLEXITY_TTCOMPLEX(self.mol_objects)
        #self['SP3CARBONS_TOTALATOM_COUNT'] = SP3CARBONS_TOTALATOM_COUNT(self.mol_objects)
        #self['SP3CARBONS_TOTALCARBON_COUNT'] = SP3CARBONS_TOTALCARBON_COUNT(self.mol_objects)
        #self['SP3CARBONS_CAR_ALLATOM_RATIO'] = SP3CARBONS_CAR_ALLATOM_RATIO(self.mol_objects)
        #self['SP3CARBONS_CHIRAL_ALLATOM_RATIO'] = SP3CARBONS_CHIRAL_ALLATOM_RATIO(self.mol_objects)
        #self['SP3CARBONS_CHIRAL_ALLCARBON_RATIO'] = SP3CARBONS_CHIRAL_ALLCARBON_RATIO(self.mol_objects)
        #self['SP3CARBONS_CHIRAL_COUNT'] = SP3CARBONS_CHIRAL_COUNT(self.mol_objects)
        #self['SP3CARBONS_CSP2_COUNT'] = SP3CARBONS_CSP2_COUNT(self.mol_objects)
        #self['SP3CARBONS_CSP2_ALLATOM_RATIO'] = SP3CARBONS_CSP2_ALLATOM_RATIO(self.mol_objects)
        #self['SP3CARBONS_CSP2_ALLCARBON_RATIO'] = SP3CARBONS_CSP2_ALLCARBON_RATIO(self.mol_objects)
        #self['SP3CARBONS_CSP3_COUNT'] = SP3CARBONS_CSP3_COUNT(self.mol_objects)
        #self['SP3CARBONS_CSP3_ALLATOM_RATIO'] = SP3CARBONS_CSP3_ALLATOM_RATIO(self.mol_objects)
        #self['SP3CARBONS_CSP3_ALLCARBON_RATIO'] = SP3CARBONS_CSP3_ALLCARBON_RATIO(self.mol_objects)
        #self['SP3CARBONS_CSP_COUNT'] = SP3CARBONS_CSP_COUNT(self.mol_objects)
        #self['SP3CARBONS_CSP_ALLATOM_RATIO'] = SP3CARBONS_CSP_ALLATOM_RATIO(self.mol_objects)
        #self['SP3CARBONS_CSP_ALLCARBON_RATIO'] = SP3CARBONS_CSP_ALLCARBON_RATIO(self.mol_objects)
        #self['RINGINFO_NUM_ALI_CARBOCYCLE'] = RINGINFO_NUM_ALI_CARBOCYCLE(self.mol_objects)
        #self['RINGINFO_NUM_ALI_HETEROCYCLE'] = RINGINFO_NUM_ALI_HETEROCYCLE(self.mol_objects)
        #self['RINGINFO_NUM_ALI_RINGS'] = RINGINFO_NUM_ALI_RINGS(self.mol_objects)
        #self['RINGINFO_NUM_ARO_CARBOCYCLE'] = RINGINFO_NUM_ARO_CARBOCYCLE(self.mol_objects)
        #self['RINGINFO_NUM_ARO_HETEROCYCLE'] = RINGINFO_NUM_ARO_HETEROCYCLE(self.mol_objects)
        #self['RINGINFO_NUM_ARO_RINGS'] = RINGINFO_NUM_ARO_RINGS(self.mol_objects)
        #self['RINGINFO_NUM_BRIDGE_ATOMS'] = RINGINFO_NUM_BRIDGE_ATOMS(self.mol_objects)
        #self['RINGINFO_NUM_SPIRO_ATOMS'] = RINGINFO_NUM_SPIRO_ATOMS(self.mol_objects)
        #self['WIENER_INDEX'] = WIENER_INDEX(self.mol_objects)
        #self['SMILES_3_2'] = SMILES_3_2(mol_smiles)
        #self['PUBCHEM_XLOGP'] = PUBCHEM_XLOGP(self.mol_objects)
        #self['PUBCHEM_TPSA'] = PUBCHEM_TPSA(self.mol_objects)
        #self['PUBCHEM_H_BOND_DONOR_COUNT'] = PUBCHEM_H_BOND_DONOR_COUNT(self.mol_objects)
        #self['PUBCHEM_H_BOND_ACCEPTOR_COUNT'] = PUBCHEM_H_BOND_ACCEPTOR_COUNT(self.mol_objects)
        #self['PUBCHEM_ROTATABLE_BOND_COUNT'] = PUBCHEM_ROTATABLE_BOND_COUNT(self.mol_objects)
        #self['PUBCHEM_HEAVY_ATOM_COUNT'] = PUBCHEM_HEAVY_ATOM_COUNT(self.mol_objects)
        #self['PUBCHEM_ATOM_STEREO_COUNT'] = PUBCHEM_ATOM_STEREO_COUNT(self.mol_objects)
        #self['PUBCHEM_DEFINED_ATOM_STEREO_COUNT'] = PUBCHEM_DEFINED_ATOM_STEREO_COUNT(self.mol_objects)
        #self['PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT'] = PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(self.mol_objects)
        #self['PUBCHEM_BOND_STEREO_COUNT'] = PUBCHEM_BOND_STEREO_COUNT(self.mol_objects)
        #self['PUBCHEM_DEFINED_BOND_STEREO_COUNT'] = PUBCHEM_DEFINED_BOND_STEREO_COUNT(self.mol_objects)
        #self['PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT'] = PUBCHEM_UNDEFINED_ATOM_STEREO_COUNT(self.mol_objects)
        #self['PUBCHEM_COVALENT_UNIT_COUNT'] = PUBCHEM_COVALENT_UNIT_COUNT(self.mol_objects)
        #self['KAPPA_SHAPE_INDEX1'] = KAPPA_SHAPE_INDEX1(self.mol_objects)
        #self['KAPPA_SHAPE_INDEX2'] = KAPPA_SHAPE_INDEX2(self.mol_objects)
        #self['KAPPA_SHAPE_INDEX3'] = KAPPA_SHAPE_INDEX3(self.mol_objects)
        #self['MCGOWAN_VOLUME'] = MCGOWAN_VOLUME(self.mol_objects)
        #self['MOE_TYPE_Labute_ASA'] = MOE_TYPE_Labute_ASA(self.mol_objects)
        #self['MOE_TYPE_PEOE_VSA'] = MOE_TYPE_PEOE_VSA(self.mol_objects)
        #self['MOE_TYPE_SMR_VSA'] = MOE_TYPE_SMR_VSA(self.mol_objects)
        #self['MOE_TYPE_SLOGP_VSA'] = MOE_TYPE_SLOGP_VSA(self.mol_objects)
        #self['MOE_TYPE_ESTATE_VSA'] = MOE_TYPE_ESTATE_VSA(self.mol_objects)
        #self['VDW_VOLUME_ABC'] = VDW_VOLUME_ABC(self.mol_objects)
        #self['ZAGREB_INDEX'] = ZAGREB_INDEX(self.mol_objects)
                                                                               #
    