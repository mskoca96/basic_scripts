#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import rdkit.Chem.rdMolAlign
from rdkit.Chem import *
from rdkit import Chem
import os




def alignTwoMol(mol1, mol2):
    #return rdkit.Chem.rdMolAlign.AlignMol(mol1, mol2)
    return rdkit.Chem.rdMolAlign.GetBestRMS(mol1, mol2)

def alignTwoMol_2(mol1, mol2):
    return rdkit.Chem.rdMolAlign.AlignMol(mol1, mol2)
    #return rdkit.Chem.rdMolAlign.GetBestRMS(mol1, mol2)

def rmsd_calc(conf, ref):
    SDFFile = conf
    suppl = SDMolSupplier(SDFFile)


    SDF_ref = ref
    suppl_ref = SDMolSupplier(SDF_ref)


    rmsd = alignTwoMol(suppl_ref[0], suppl[0])
    return rmsd


def rmsd_calc_2(conf, ref):
    SDFFile = conf
    suppl = SDMolSupplier(SDFFile)


    SDF_ref = ref
    suppl_ref = SDMolSupplier(SDF_ref)


    rmsd = alignTwoMol_2(suppl_ref[0], suppl[0])
    return rmsd

def rmsd_matrix(struct_dir):
    from scipy.cluster.hierarchy import ward, fcluster
    from scipy.spatial.distance import pdist
    all_files = os.listdir(struct_dir)
    d = pd.DataFrame(0, index=all_files, columns=all_files).astype(float)
    for i in range(0,len(all_files)):
        for j in range(0,len(all_files)):
            mol1 = struct_dir + "/"+ all_files[i]
            mol2 = struct_dir + "/"+ all_files[j]
            rmsd = rmsd_calc(mol1, mol2)
            d.iloc[i,j] = rmsd
            print(rmsd)
    Z = ward(pdist(d))
    clustering = (fcluster(Z, t=0.5, criterion='distance'))
    d["cluster"] = clustering
    d.to_csv("Results_%s.csv"%(str(struct_dir)))




import argparse

def mainAlingTwoMol():
    parser = argparse.ArgumentParser(description="Give something ...")
    parser.add_argument("-i", "--mol_path", type=str, required=True, help="")
    parser.add_argument("-r", "--mol_ref_path", type=str, required=True, help="")
    parser.add_argument("-o", "--output", type=str, required=True, help="")
    args=parser.parse_args()
    mol_path=args.mol_path
    mol_ref_path=args.mol_ref_path
    output=args.output
    
    SDFFile = mol_path
    suppl = SDMolSupplier(SDFFile, removeHs= False)

    
    suppl_h = SDMolSupplier(SDFFile) #removeHs= False)

    SDF_ref = mol_ref_path
    suppl_ref = SDMolSupplier(SDF_ref, removeHs=False)


    rmsd = alignTwoMol_2(suppl_ref[0], suppl[0])
    #rmsd = rmsd_calc_2(mol_path, mol_ref_path)
    print("RMSD-->", rmsd)
    writer = rdkit.Chem.rdmolfiles.SDWriter(output)
    writer.write(suppl_h[0])

mainAlingTwoMol()
#rmsd_matrix("deneme")
