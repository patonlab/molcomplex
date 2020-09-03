# rdkit
import rdkit
from rdkit import Chem

#import the descriptors functions

class mol:
    #definition of molobjects
    def __init__(self,smi):

        #general names for molecule
        self.input_smiles = smi
        self.mol = Chem.MolFromSmiles(smi)
        self.canonical_smiles = Chem.rdmolfiles.MolToSmiles(self.mol)
        self.inchi = Chem.inchi.MolToInchi(self.mol)
        self.inchikey =  Chem.inchi.MolToInchiKey(self.mol)
        # self.iupac_name

        #molecular pro[erties
        self.molecular_formula = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
        self.molecular_weight = Chem.Descriptors.ExactMolWt(self.mol)

        #assing functions for mol Descriptors from descriptors .py

if __name__ == "__main__":

    #just for testing the class is working
    smi = 'CCCC'
    mol_object = mol(smi)
    print('SMILES',mol_object.smiles)
    print('MOLFromSMILES',mol_object.mol)
