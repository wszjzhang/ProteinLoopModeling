"""Combine template and input structure with procutes analysis

   This module reads in the loop structure template
   fit the loop into the incomplete PDB and write complete PDB

   version 0.1 only consider loop in the middle of the structure
"""
    
import sys
import os
import shutil


import numpy as np


from procrustes import procrustes
from getTemplate import searchDB
from getTemplate import getTemp


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
    inpdb = np.array([line.split() for line in fil if 'CA' in line])
    fil.close()

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
        lpcrd = np.array([line.split() for line in fil][1:-1])
        fil.close()
        lpcrds.append(lpcrd[:,6:9].astype(np.float))

    pdbStart = resids[0]
    pdbEnd   = resids[-1] 
    aalength = len(seq)
    
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
    
    #fit with 3 amino acids of each end
    iOverlapStart1 = int(iLoopStart) - 3 
    iOverlapStart2 = int(iLoopStart)  
    iOverlapEnd1 = int(iLoopEnd) + 1 
    iOverlapEnd2 = int(iLoopEnd) + 4 

    #transform lpcrd, fit in incomplete structure
    strCompare = np.concatenate((incrd[iOverlapStart1:iOverlapStart2],
        incrd[iOverlapEnd1:iOverlapEnd2]))
    lpCompare = np.concatenate((lpcrd[:3], lpcrd[-3:]))
    score,xyz_tf,tf = procrustes(strCompare, lpCompare,False)
    lpcrd_recon = np.dot(lpcrd, tf['rotation']) + np.tile(tf['translation'][0],[lpcrd.shape[0],1]);

    #combine lpcrd_recon with inpdb
    completeStructure = np.concatenate((incrd[:iLoopStart], lpcrd_recon[4:-3], incrd[iLoopStart:]))
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


def main():
    # Get the name from the command line
    if len(sys.argv) == 3:
        #log = sys.argv[1]
        print 'reading input PDB and sequence file'
        print '......\n'
    else:
        print 'Usage: ./loopmodeling.py seq.fsa target.pdb'
        return
    
    pdbfile = sys.argv[1]
    seqfile = sys.argv[2]

    loopfiles = getLoopfiles()
    target,seq,incrd,lpcrds,pdbStart,loopStart,loopEnd = readParameters(seqfile,pdbfile,loopfiles)
    completecrd = combStructure(incrd,lpcrd,loopStart,loopEnd)
    writePDB(target, seq, pdbStart, completecrds)
    
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
        main()
