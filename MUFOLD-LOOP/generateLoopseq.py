#!/usr/bin/env python -tt
# Copyright 2015 Jiong Zhang

"""read in the full seq file and input structure pdb
   write loop seqfile = 3AA+loop+3AA
"""


def genloopseq(seqfile, pdbfile):

    #get reesidue ids in pdbfile
    fil = open(pdbfile,'r')
    resids = [line[22:27] for line in fil if 'CA' in line]
    fil.close()
    resids = map(int, resids)
    resids.sort()

    #full sequence file: seqfile
    fil = open(seqfile,'r')
    seq = [line for line in fil][-1]
    if not seq[-1].isalpha():
        seq = seq[:-1]
    fil.close()

    if resids[-1]-resids[0]+1 != len(seq):
        print "sequence file and pdb file do not match !!!"
        return

    #print loop seq file
    loopres = list(set(range(resids[0],resids[-1]+1)) - set(resids))
    loopres.sort()
    StartId = loopres[0]-resids[0]-3
    print loopres[0]
    print resids[0]
    lplen = len(loopres)
    fil = open('./loopseq.fsa','w')
    fil.write('> 3AA+loop+3AA \n')
    fil.write(seq[StartId:StartId+lplen+6]+'\n')
    fil.close()

    
