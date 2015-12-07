#!/usr/bin/env python -tt
#copyright 2015 Jiong Zhang

"""This file is the main script to call FastLoop 
   fit the loop into the incomplete PDB and write complete PDB
 
   version 0.1 only consider loop in the middle of the structure
"""

print "FastLoop start"

import sys

from getTemplate import searchDB
from getTemplate import getTemp

from combineStructure import getLoopfiles
from combineStructure import readParameters

from generateModel import generate_model


def main():
   # Get the name from the command line
    if len(sys.argv) == 4 or len(sys.argv) == 5:
    #log = sys.argv[1]
        print 'reading input PDB and sequence file'
        print '......\n'
    else:
        print 'Usage: ./FastLoop.py fullseq.fsa loopseq.fsa target.pdb [savesteps]'
        return

    fullseqfile = sys.argv[1]
    loopseqfile = sys.argv[2]
    pdbfile = sys.argv[3]
    
    if len(sys.argv) == 5:
        savesteps = sys.argv[4]
    else:
        savesteps = False

    loop_length = searchDB(loopseqfile)
    
    fil = open(loopseqfile,'r')
    loopseq = [line for line in fil][-1]
    if not loopseq[-1].isalpha():
        loopseq = loopseq[:-1]
    fil.close()
    print loopseq

    getTemp("blastout.txt", loop_length)

    loopfiles = getLoopfiles()
    #loopfiles = ["./templates/3E59_loop5.pdb"]

    target, seq, incrd, lpcrds, pdbStart, loopStart, loopEnd = readParameters(fullseqfile, pdbfile, loopfiles)
    print "==============================================================="
    print "Modeling loop for",target,"from residue",loopStart,"to",loopEnd 
    print "==============================================================="

    for loopn in range(len(lpcrds)):
        generate_model(pdbfile, incrd, lpcrds[loopn], loopseq, loopStart, loopEnd, target, loopn, savesteps)

    # QA 
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
