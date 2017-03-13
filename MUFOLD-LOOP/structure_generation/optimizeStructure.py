#!/usr/bin/env python -tt
#copyright 2015 Jiong Zhang

"""This module read in a structure
   1 generate distance matrix of the structure
   2 assign all value along the diagnal to 3.8
   3 use shortest path to reassign pair distance values
   4 use MDS convert distance matrix to 3D model
   5 combine loop part with initial incomplete PDB
   6 repeat 1-5 for N times
"""


import numpy as np
from scipy.spatial import distance


from Dijkstra import Dijkstra
from classicalMDS import cmdscale
from sklearn import manifold
from procrustes import procrustes

def pdb2DM(pdbfile):

    print "processing: ", pdbfile
    fil = open(pdbfile,'r')
    inpdb = np.array([line.split() for line in fil if 'CA' in line])
    fil.close()

#    resids = map(int,list(inpdb[:,5]))

    incrd = inpdb[:,5:8].astype(np.float)
    dm = distance.squareform(distance.pdist(incrd,"euclidean"))
    return dm


def crd2DM(incrd):
    dm = distance.squareform(distance.pdist(incrd,"euclidean"))
    return dm


def shortest_path(dm, lpStart, lpEnd):
    success = True
    # reassign adjacent atom pair distance bigger than 3.9 to 3.8
    for i in range(lpStart,lpEnd+2):
        if dm[i-1][i] > 3.9:
            print "pair distance too big:", i, i-1, dm[i-1][i] 
            dm[i-1][i] = 3.8
            dm[i][i-1] = 3.8
            success = False
        if dm[i-1][i] < 3.7:
            print "pair distance too small:", i, i-1, dm[i-1][i] 
            dm[i-1][i] = 3.8
            dm[i][i-1] = 3.8
            success = False

    # reassign pair distances with shortest path
    nCA =len (dm)
    Dijk = False
    if Dijk == False:
        n = 1  # number of dm optimization rounds
        for ni in range(n):
            for k in range(1,nCA):          #    k: third residue
                for i in range(lpStart,lpEnd+1):      #    i: residue in loop region 
                    for j in range(0,i-1)+range(i+2,nCA):#    j: distant residue to i  
                        if dm[i,j] > dm[i,k]+dm[k,j]:
                            dm[i,j] = dm[i,k]+dm[k,j];
                            dm[j,i] = dm[i,j]
#                           print i, j
        return dm,success
    else:
        for i in range(len(dm)):
            for j in range(i+1, len(dm)):
                dm[i][j] = Dijkstra(dm, i, j)
                dm[j][i] = dm[i][j]
        return dm,success


def MDSdm2crd(incrd, dm):
    sk = True
    if sk:
        # MDS distanceM => coordinates
        mds = manifold.MDS(n_components=3, max_iter=3000, eps=1e-9,
                   dissimilarity="precomputed", n_jobs=1)
        new_pos = mds.fit(dm).embedding_
        d,crd,tform = procrustes(incrd, new_pos, False)
    else:
        new_pos = cmdscale(dm)[0][:,0:3]
        d,crd,tform = procrustes(incrd, new_pos, False)
    return crd



