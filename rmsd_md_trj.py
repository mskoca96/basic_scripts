#!/usr/bin/env python
# -*- coding: utf-8 -*-


from math import *

def rmsd_to_atom(pdb_file, atom_1, atom_2):
    trj_atom_1 = []
    trj_atom_2 = []
    pdb = open(pdb_file,"r")
    for line in pdb:
        try :
            if int(line[6:11].strip()) == atom_1:
                trj_atom_1.append(line[31:56])
            if int(line[6:11].strip()) == atom_2:
                trj_atom_2.append(line[31:56])
        except:
            pass
    output = open("rmsd_between_%d_%d.csv"%(atom_1, atom_2),"w")
    for i in range(0,len(trj_atom_1)):
        coord_1 = trj_atom_1[i].split()
        coord_2 = trj_atom_2[i].split()
        x = pow((float(coord_1[0])-float(coord_2[0])),2)
        y = pow((float(coord_1[1])-float(coord_2[1])),2)
        z = pow((float(coord_1[2])-float(coord_2[2])),2)
        rmsd = sqrt(x+y+z)
        output.write((str(i+1)+","+str(rmsd)+"\n"))



rmsd_to_atom("deneme.pdb",1,2)

