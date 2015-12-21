#!/usr/bin/env python -tt
#copyright 2015 Jiong Zhang

"""Use combineStructure.py and optimizeStructure.py generate final model
   
   generate distance matrix of the structure
   assign all value along the diagnal to 3.8
   use shortest path to reassign pair distance values
   use MDS convert distance matrix to 3D model
   combine loop part with initial incomplete PDB
   repeat 1-5 for N times
"""

import os
import numpy as np

from procrustes import procrustes

from combineStructure import combStructure
from combineStructure import convertAA

from optimizeStructure import crd2DM
from optimizeStructure import shortest_path
from optimizeStructure import MDSdm2crd


def generate_model(pdbfile, incrd, lpcrd, loopseq, loopStart, loopEnd, target, loopn, savesteps):
    # combine loop with input structure and optimize complete structure
    optimization_done = False
    r = 1
    while not optimization_done:
        # iterate until all adjascent Ca-Ca dist around 3.8 
        print "Optimization round ", r
        r += 1
        # combine loop structure with input pdb coordinates
        completecrds = combStructure(incrd, lpcrd, loopStart, loopEnd)
        # convert coordinates to distance matrix
        dist_matrix = crd2DM(completecrds)
        # optimize structure with shortest path and MDS
        dist_matrix,optimization_done = shortest_path(dist_matrix)
        optimized_crd = MDSdm2crd(completecrds, dist_matrix)

        # update loop structure with optimized structure
        lpcrd = optimized_crd[loopStart-3:loopEnd+4]

        if savesteps:
            np.savetxt("intermediate_"+str(r)+".crd",optimized_crd)
        if r > 200:
            optimization_done = True
    # put lpcrd in right orientation
    #transform lpcrd, fit in incomplete structure
    strCompare = incrd[loopStart-3:loopStart+3]             #anchor res: ls-3, ls-2, ls-1, ls,ls+1,ls+2
    lpCompare = np.concatenate((lpcrd[:3], lpcrd[-3:]))
    score,xyz_tf,tf = procrustes( lpCompare, strCompare, False)
    lpcrd_recon = np.mat(lpcrd)*np.mat(tf['rotation'])+np.mat(tf['translation']);
    lpcrd_recon = np.array(lpcrd_recon)

    ### write optimized loop pdb

    # read in initial full atom model
    fil = open(pdbfile,'r')
    inpdb = np.array([[line[0:4],line[4:11],line[13:17],line[17:20],line[20:22],line[22:27], line[30:38],line[38:46],line[46:54]] for line in fil if 'ATOM' in line])
    fil.close()
    resids = map(int,list(inpdb[:,5]))
    
    # get loop resid
    pdbStart = resids[0]
    pdbEnd   = resids[-1]
    loopResids = sorted(list(set(range(pdbStart,pdbEnd+1))-set(resids)))
    loopSeq = loopseq[3:-3]

    # obtain the full atom structure of the loop
    loopname = target+'loop'+str(loopn)+'.pdb'
    print "writing ", loopname, "..."
    fil = open(loopname,'w')
    print loopSeq, len(loopSeq), len(loopResids), len(lpcrd_recon[3:-3])
    # write loop structure
    for k in range(len(loopResids)):
        fil.write('ATOM  %5d  CA  %3s %5d    %8.3f%8.3f%8.3f  1.00  0.50\n'
                    %(k, convertAA(loopSeq[k]), int(loopResids[k]), float(lpcrd_recon[k+3][0]), float(lpcrd_recon[k+3][1]), float(lpcrd_recon[k+3][2])))
    fil.close()
    # call PULCHRA add side chain
    PULCHRA = '/Users/jiong/jiong/projects/5_MUFOLD-FL/MUFOLD-FL_py_skMDS/pulchra304/bin/osx/pulchra '
    pulchracmd = PULCHRA + loopname
    os.system(pulchracmd)
    

    ####################################################################
    ###  combine fullatom input structure and fullatom loo structure ###
    ####################################################################


    # read in full atom loop structure
    fil = open(target+'loop'+str(loopn)+'.rebuilt.pdb','r')
    looppdb = np.array([[line[0:4],line[4:11],line[13:17],line[17:20],line[20:22],line[22:27],line[30:38],line[38:46],line[46:54]] for line in fil if 'ATOM' in line])
    fil.close()
      
    # combine loop and input structure: write complete structure
    name = './complete_strs/'+target+'complete'+str(loopn)+'.pdb'
    print "writing ", name, "..."
    fil = open(name,'w')
    # structure before loop
    i = 0
    while int(inpdb[i][5]) < loopResids[0]:
        fil.write('ATOM  %5d  %4s%3s %5d    %8.3f%8.3f%8.3f  1.00  0.50\n'
                    %(i+1, inpdb[i][2], inpdb[i][3], int(inpdb[i][5]), float(inpdb[i][6]), float(inpdb[i][7]), float(inpdb[i][8])))
        i += 1
    j = i # track input pdb index 
    # loop structure
    for k in range(len(looppdb)):
        fil.write('ATOM  %5d  %4s%3s %5d    %8.3f%8.3f%8.3f  1.00  0.50\n'
                    %(i+1, looppdb[k][2], looppdb[k][3], int(looppdb[k][5]), float(looppdb[k][6]), float(looppdb[k][7]), float(looppdb[k][8])))
        i += 1
    # structure after loop
    while len(inpdb[j:]) > 0:
        fil.write('ATOM  %5d  %4s%3s %5d    %8.3f%8.3f%8.3f  1.00  0.50\n'
                   %(i+1, inpdb[j][2], inpdb[j][3], int(inpdb[j][5]), float(inpdb[j][6]), float(inpdb[j][7]), float(inpdb[j][8])))
        j += 1
        i += 1

    fil.close()

    # record Nround
    fil = open("optimization.log",'a')
    fil.write("%s optimized in rounds: %d \n" % (name, r))
    fil.close()


    # clean intermediate files
    os.remove(target+'loop'+str(loopn)+'.pdb')
    os.remove(target+'loop'+str(loopn)+'.rebuilt.pdb')

