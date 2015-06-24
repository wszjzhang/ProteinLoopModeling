#!/usr/bin/python -tt
#copyright 2015 Jiong Zhang

"""This module read in a structure
   1 generate distance matrix of the structure
   2 assign all value along the diagnal to 3.8
   3 use shortest path to reassign pair distance values
   4 use MDS convert distance matrix to 3D model
   5 combine loop part with initial incomplete PDB
   6 repeat 1-5 for 3 times
"""


import numpy as np
import scipy


from scipy.spatial import distance


from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA


from Dijkstra import Dijkstra
from classicalMDS import cmdscale

def pdb2DM(pdbfile):

    print "processing: ", pdbfile
    fil = open(pdbfile,'r')
    inpdb = np.array([line.split() for line in fil if 'CA' in line])
    fil.close()

    resids = map(int,list(inpdb[:,5]))

    incrd = inpdb[:,5:8].astype(np.float)
    print "calculate distance matrix"
    dm = distance.squareform(distance.pdist(incrd,"euclidean"))
    return dm


def crd2DM(incrd):
    print "calculating distance matrix"
    dm = distance.squareform(distance.pdist(incrd,"euclidean"))
    return dm


def shortest_path(dm, Dijkstra):
    nCA = len(dm)
    # reassign adjacent atom pair distance bigger than 3.9 to 3.8
    for i in range(1,nCA):
        if dm[i-1][i] > 3.9:
            dm[i-1][i] = 3.8
            dm[i][i-1] = 3.8
            print i
    # reassign pair distances with shortest path
    if Dijkstra == False:
        for k in range(1,nCA):
            for i in range(1,nCA):
                for j in range(i+1,nCA):
                    if dm[i,j]>dm[i,k]+dm[k,j]:
                        dm[i,j]=dm[i,k]+dm[k,j];
                        dm[j,i]= dm[i,j]
#                        print i, j
        return dm
    else:
        for i in range(len(dm)):
            for j in range(i+1, len(dm)):
                dm[i][j] = Dijkstra(G, i, j)
                dm[j][i] = dm[i][j]
        return dm


def MDSdm2crd(dm):
    crd = cmdscale(dm)[0][:,0:3]
    return crd



