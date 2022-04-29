#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import rdkit.Chem.rdMolAlign
from rdkit.Chem import *
from rdkit import Chem
import os


def alignTwoMol_2(mol1, mol2):
    return rdkit.Chem.rdMolAlign.CalcRMS(mol1, mol2)
    #return rdkit.Chem.rdMolAlign.GetBestRMS(mol1, mol2)


def rmsd_calc(conf, ref):
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
    clustering = (fcluster(Z, t=30, criterion='distance'))
    d["cluster"] = clustering
    d.to_csv("Results_%s.csv"%(str(struct_dir)))

rmsd_matrix("ligands")
