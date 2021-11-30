#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	    	  used for genrating a proudfoot index    	 #
#####################################################.

import rdkit
from rdkit import Chem
import numpy as np


def calc_ca(num_uni_detail_list_1,num_uni_detail_list_2,total_all_details):
    if total_all_details == 0: return 0
    ca = 0
    for list_1 in num_uni_detail_list_1:
        ca -= (list_1[1]/total_all_details)*(np.log2(list_1[1]/total_all_details))
    for list_2 in num_uni_detail_list_2:
        ca -= (list_2[1]/total_all_details)*(np.log2(list_2[1]/total_all_details))

    ca += np.log2(total_all_details)
    return ca

def calc_cm_star(ca_all):
    sum = 0
    for ca in ca_all:
        sum += 2**ca
    return np.log2(sum)

def proudfoot_index(mol):

    ca_all = []
    for root_atom in range(0,mol.GetNumAtoms()):
        all_detail_list_1, all_detail_list_2 = [],[]
        len_num_bonds_1,len_num_bonds_2  =  None,None

        type_list_1, type_list_2 = [], []
        if mol.GetAtomWithIdx(root_atom).GetSymbol() != 'H':
            #getting all paths of lengrth 1 and 2
            num_bonds_1 = Chem.rdmolops.FindAllPathsOfLengthN(mol,1,rootedAtAtom=root_atom,useHs=True)
            num_bonds_2 = Chem.rdmolops.FindAllPathsOfLengthN(mol,2,rootedAtAtom=root_atom,useHs=True)

            len_num_bonds_1 = len(list(num_bonds_1))

            #only taking Hydrogen for path length 1 and making a list of bond indexes for path length 1 and 2
            final_1_bond, final_2_bond = [],[]
            for i,bond_1 in enumerate(reversed(num_bonds_1)):

                bond = list(bond_1)[0]
                atom_1 = mol.GetBondWithIdx(bond).GetBeginAtomIdx()
                atom_2 = mol.GetBondWithIdx(bond).GetEndAtomIdx()
                if atom_1 == root_atom:
                    atom_remove = atom_2
                if atom_2 == root_atom:
                    atom_remove = atom_1

                if mol.GetAtomWithIdx(atom_remove).GetSymbol() != 'H':
                    len_num_bonds_1 -= 1
                else:
                    final_1_bond.append(bond)

            for i,bond_2 in enumerate(reversed(num_bonds_2)):
                final_2_bond.append(list(bond_2))

            len_num_bonds_2 = len(list(num_bonds_2))

            total_for_each_atom = len_num_bonds_1 + len_num_bonds_2

            #assigning the atom descriptions for each bond bond distance 1
            for each_1 in final_1_bond:
                atom_1_no = mol.GetBondWithIdx(bond).GetBeginAtomIdx()
                atom_2_no = mol.GetBondWithIdx(bond).GetEndAtomIdx()

                atom_1_exp_non_H_connection = 0
                #ATOM_1 DETAILS
                atom_1 = mol.GetAtomWithIdx(atom_1_no)
                atom_1_atom_num = atom_1.GetAtomicNum()
                atom_1_tot_connection = atom_1.GetTotalValence()
                neighbours_1 = atom_1.GetNeighbors()
                for niegh in neighbours_1:
                    if niegh.GetSymbol() != 'H':
                        atom_1_exp_non_H_connection +=1

                atom_2_exp_non_H_connection = 0
                #ATOM_1 DETAILS
                atom_2 = mol.GetAtomWithIdx(atom_2_no)
                atom_2_atom_num = atom_2.GetAtomicNum()
                atom_2_tot_connection = atom_2.GetTotalValence()
                neighbours_2 = atom_2.GetNeighbors()
                for niegh in neighbours_2:
                    if niegh.GetSymbol() != 'H':
                        atom_2_exp_non_H_connection +=1

                all_detail_list_1.append([atom_1_atom_num,atom_1_tot_connection,atom_1_exp_non_H_connection, atom_2_atom_num,atom_2_tot_connection,atom_2_exp_non_H_connection])


            #assigning the atom descriptions for each bond distance 2
            for each_2 in final_2_bond:
                atom_1_no = mol.GetBondWithIdx(each_2[0]).GetBeginAtomIdx()
                atom_2_no = mol.GetBondWithIdx(each_2[0]).GetEndAtomIdx()

                if atom_2_no == root_atom:
                    temp = atom_1_no
                    atom_1_no = root_atom
                    atom_2_no = temp

                assign_3_no = mol.GetBondWithIdx(each_2[1]).GetBeginAtomIdx()
                assign_4_no = mol.GetBondWithIdx(each_2[1]).GetEndAtomIdx()

                if assign_3_no == atom_1_no or assign_3_no == atom_2_no :
                    atom_3_no = assign_4_no

                if assign_4_no == atom_1_no or assign_4_no == atom_2_no :
                    atom_3_no = assign_3_no

                atom_1_exp_non_H_connection = 0
                #ATOM_1 DETAILS
                atom_1 = mol.GetAtomWithIdx(atom_1_no)
                atom_1_atom_num = atom_1.GetAtomicNum()
                atom_1_tot_connection = atom_1.GetTotalValence()
                neighbours_1 = atom_1.GetNeighbors()
                for niegh in neighbours_1:
                    if niegh.GetSymbol() != 'H':
                        atom_1_exp_non_H_connection +=1

                atom_2_exp_non_H_connection = 0
                #ATOM_2 DETAILS
                atom_2 = mol.GetAtomWithIdx(atom_2_no)
                atom_2_atom_num = atom_2.GetAtomicNum()
                atom_2_tot_connection = atom_2.GetTotalValence()
                neighbours_2 = atom_2.GetNeighbors()
                for niegh in neighbours_2:
                    if niegh.GetSymbol() != 'H':
                        atom_2_exp_non_H_connection +=1

                atom_3_exp_non_H_connection = 0
                #ATOM_3 DETAILS
                atom_3 = mol.GetAtomWithIdx(atom_3_no)
                atom_3_atom_num = atom_3.GetAtomicNum()
                atom_3_tot_connection = atom_3.GetTotalValence()
                neighbours_3 = atom_3.GetNeighbors()
                for niegh in neighbours_3:
                    if niegh.GetSymbol() != 'H':
                        atom_3_exp_non_H_connection +=1

                all_detail_list_2.append([atom_1_atom_num,atom_1_tot_connection,atom_1_exp_non_H_connection, atom_2_atom_num,atom_2_tot_connection,atom_2_exp_non_H_connection, atom_3_atom_num,atom_3_tot_connection,atom_3_exp_non_H_connection])


            total_all_details = len(all_detail_list_1) + len(all_detail_list_2)

            uni_detail_list_1,num_uni_detail_list_1,uni_detail_list_2,num_uni_detail_list_2 = [],[],[],[]
            #unique of atom 1
            for list_1 in all_detail_list_1:
                if list_1 in uni_detail_list_1:
                    for uni_list_1 in num_uni_detail_list_1:
                        if uni_list_1[0] == list_1:
    #                         print(uni_list_1)
                            uni_list_1[1] += 1
                else:
                    uni_detail_list_1.append(list_1)
                    num_uni_detail_list_1.append([list_1,1])

            #unique of atom 1
            for list_2 in all_detail_list_2:
                if list_2 in uni_detail_list_2:
                    for uni_list_2 in num_uni_detail_list_2:
                        if uni_list_2[0] == list_2:
                            uni_list_2[1] += 1
                else:
                    uni_detail_list_2.append(list_2)
                    num_uni_detail_list_2.append([list_2,1])

            ca = calc_ca( num_uni_detail_list_1, num_uni_detail_list_2, total_all_details)

            ca_all.append(ca)

    cm = np.sum(ca_all)

    cm_star = calc_cm_star(ca_all)

    return cm,cm_star
