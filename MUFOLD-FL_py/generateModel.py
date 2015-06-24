#!/usr/bin/python -tt
#copyright 2015 Jiong Zhang

"""Use combineStructure.py and optimizeStructure.py generate final model
   
   generate distance matrix of the structure
   assign all value along the diagnal to 3.8
   use shortest path to reassign pair distance values
   use MDS convert distance matrix to 3D model
   combine loop part with initial incomplete PDB
   repeat 1-5 for 3 times
"""

import os
import numpy as np


from combineStructure import getLoopfiles
from combineStructure import readParameters
from combineStructure import combStructure
from combineStructure import writePDB
from combineStructure import convertAA

from optimizeStructure import crd2DM
from optimizeStructure import pdb2DM
from optimizeStructure import shortest_path
from optimizeStructure import MDSdm2crd


def generate_model(pdbfile, incrd, lpcrd, loopseq, loopStart, loopEnd, target, loopn):
    # combine loop with input structure and optimize complete structure
    for r in range(10):
        # 3 iteration 
        print "Optimization round ", r+1
        # combine loop structure with input pdb coordinates
        completecrds = combStructure(incrd, lpcrd, loopStart, loopEnd)
        # convert coordinates to distance matrix
        dist_matrix = crd2DM(completecrds)
        # optimize structure with shortest path and MDS
        dist_matrix = shortest_path(dist_matrix, False)
        optimized_crd = MDSdm2crd(dist_matrix)
        # update loop structure with optimized structure
        lpcrd = optimized_crd[loopStart-3:loopEnd+4]
    

    # build full atom model 
    # read in initial full atom model
    fil = open(pdbfile,'r')
    inpdb = np.array([line.split() for line in fil if 'ATOM' in line])
    fil.close()

    resids = map(int,list(inpdb[:,5]))
    pdbStart = resids[0]
    pdbEnd   = resids[-1]

    # get loop resid
    loopResids = sorted(list(set(range(pdbStart,pdbEnd+1))-set(resids)))
    loopSeq = loopseq[3:-3]


    # obtain the full atom structure of the loop
    loopname = target+'loop'+str(loopn)+'.pdb'
    print "writing ", loopname, "..."
    fil = open(loopname,'w')
    print loopSeq, len(loopSeq), len(loopResids), len(lpcrd[3:-3])
    # loop structure
    for k in range(len(loopResids)):
        fil.write('ATOM  %5d  CA  %3s %5d    %8.3f%8.3f%8.3f  1.00  0.50\n'
                    %(k, convertAA(loopSeq[k]), loopResids[k], lpcrd[k+3][0]    , lpcrd[k+3][1], lpcrd[k+3][2]))
    fil.close()
    # call PULCHRA
    PULCHRA = '/Volumes/JiongData/MUFOLD-FL/MUFOLD-FL0_0/pulchra304/bin/osx/pulchra '
    pulchracmd = PULCHRA + loopname
    os.system(pulchracmd)
    
    # read in full atom loop structure
    fil = open(target+'loop'+str(loopn)+'.rebuilt.pdb','r')
    looppdb = np.array([line.split() for line in fil if 'ATOM' in line])
    fil.close()
    
    # write complete structure
    name = target+'complete'+str(loopn)+'.pdb'
    print "writing ", name, "..."
    fil = open(name,'w')
    # structure before loop
    i = 0
    while int(inpdb[i][5]) < loopResids[0]:
        fil.write('ATOM  %5d  %3s %3s %5s    %8s%8s%8s  1.00  0.50\n'
                    %(i+1, inpdb[i][2], inpdb[i][3], inpdb[i][5], inpdb[i][6], inpdb[i][7], inpdb[i][8]))
        i += 1
    j = i # track input pdb index 

    # loop structure
    for k in range(len(looppdb)):
        fil.write('ATOM  %5d  %3s %3s %5s    %8s%8s%8s  1.00  0.50\n'
                    %(i+1, looppdb[k][2], looppdb[k][3], looppdb[k][4], looppdb[k][5], looppdb[k][6], looppdb[k][7]))
        i += 1
    # structure after loop
    while len(inpdb[j:]) > 0:
        fil.write('ATOM  %5d  %3s %3s %5s    %8s%8s%8s  1.00  0.50\n'
                   %(i+1, inpdb[j][2], inpdb[j][3], inpdb[j][5], inpdb[j][6], inpdb[j][7], inpdb[j][8]))
        j += 1
        i += 1

    fil.close()

    # clean intermediate files
    os.remove(target+'loop'+str(loopn)+'.pdb')
    os.remove(target+'loop'+str(loopn)+'.rebuilt.pdb')

