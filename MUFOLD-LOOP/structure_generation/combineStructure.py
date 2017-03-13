"""Combine template and input structure with procutes analysis

   This module reads in the loop structure template
   fit the loop into the incomplete PDB and write complete PDB

   version 0.1 only consider single loop in the middle of the structure
"""
    
import os

import numpy as np

from procrustes import procrustes


def readParameters(seqfile, pdbfile, lpfiles):
    """read parameters and obtain arrays for objs
       incrd:       input coordinate,
       lpcrds:      coordinates of loop structure,
       pdbStart:    start resid in input PDB,
       loopStart:   loop start residue number,
       loopEnd:     loop end residue number ,
       seq:         full sequence AA string
    """

    target = pdbfile.split('.')[0]
    
    print 'processing target', target 
    fil = open(pdbfile,'r')
    inpdb = np.array([[line[0:4],line[4:11],line[13:17],line[17:20],line[20:22],line[22:27],line[30:38],line[38:46],line[46:54]] for line in fil if 'CA' in line])
    fil.close()
    
    #full sequence file: seqfile
    fil = open(seqfile,'r')
    seq = [line for line in fil][-1]
    if not seq[-1].isalpha():
        seq = seq[:-1]
    fil.close()
    
    resids = map(int,list(inpdb[:,5]))

    incrd = inpdb[:,6:9].astype(np.float)

    lpcrds = []
    for lpfile in lpfiles:
        fil = open(lpfile,'r')
        lpcrd = np.array([[line[0:4],line[4:11],line[13:17],line[17:20],line[20:22],line[22:27],line[30:38],line[38:46],line[46:54]] for line in fil if 'CA' in line])
        fil.close()
        lpcrds.append(lpcrd[:,6:9].astype(np.float))

    pdbStart = resids[0]
    pdbEnd   = resids[-1] 
    
    # loopStart 0 base from first residue in input pdb
    loopResids = sorted(list(set(range(pdbStart,pdbEnd+1))-set(resids)))
    loopStart = loopResids[0] - pdbStart
    loopEnd = loopResids[-1] - pdbStart
    return target,seq,incrd,lpcrds,pdbStart,loopStart,loopEnd


def getLoopfiles():
    loopfiles = []
    for dirname , dirnames, filenames in os.walk('./templates'):
        for filename in filenames:
            loopfiles.append(os.path.join(dirname, filename))
    return loopfiles


def combStructure(incrd,lpcrd,iLoopStart,iLoopEnd):
    """this function combine loop coordinates 
        with incomplete PDB coordinates
        and obtain complete PDB
    """
    
    #in incrd loop nodt included, thus overlap parts are continuous in incrd
    #transform lpcrd, fit in incomplete structure
    strCompare = incrd[iLoopStart-3:iLoopStart+3]
    lpCompare = np.concatenate((lpcrd[:3], lpcrd[-3:]))
    score,xyz_tf,tf = procrustes(strCompare, lpCompare, False, False)
    lpcrd_recon = np.mat(lpcrd)*np.mat(tf['rotation'])+np.mat(tf['translation']);
    lpcrd_recon = np.array(lpcrd_recon)
    #combine lpcrd_recon with inpdb
    completeStructure = np.concatenate((incrd[:iLoopStart], lpcrd_recon[3:-3], incrd[iLoopStart:]))
    return completeStructure


def convertAA(s):
    codeAA = {'A':'ALA',
              'R':'ARG',
              'N':'ASN',
              'D':'ASP',
              'B':'ASX',
              'C':'CYS',
              'E':'GLU',
              'Q':'GLN',
              'Z':'GLX',
              'G':'GLY',
              'H':'HIS',
              'I':'ILE',
              'L':'LEU',
              'K':'LYS',
              'M':'MET',
              'F':'PHE',
              'P':'PRO',
              'S':'SER',
              'T':'THR',
              'W':'TRP',
              'Y':'TYR',
              'V':'VAL'}
    return codeAA[s]


def writePDB(target, seq, pdbStart, structures):
    """get sequence string and coordinates
        write pdb files
    """
    loopn = 1
    for structure in structures:
        name = target+'complete'+str(loopn)+'.pdb'
        resid = pdbStart
        #if os.path.isfile(name):
        #    os.remove(name)
        print "writing ", name, "..."
        fil = open(name,'w')
        for k in xrange(len(structure)):
            fil.write('ATOM  %5d  CA  %3s %5d    %8.3f%8.3f%8.3f  1.00  0.50\n' 
                   %(resid, convertAA[seq[k]], resid, structure[k][0], structure[k][1], structure[k][2]))
            resid += 1
        fil.close()
        
        loopn += 1


