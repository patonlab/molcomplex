# rdkit
import rdkit
from rdkit import Chem

#import the descriptors functions

class mol:
    #definition of molobjects
    def __init__(self,smi):
        self.smiles = smi
        self.mol = Chem.MolFromSmiles(smi)

        #assing functions for mol Descriptors from descriptors .py

if __name__ == "__main__":

    #just for testing the class is working
    smi = 'CCCC'
    mol_object = mol(smi)
    print('SMILES',mol_object.smiles)
    print('MOLFromSMILES',mol_object.mol)
