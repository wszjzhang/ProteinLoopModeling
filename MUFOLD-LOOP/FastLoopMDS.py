#!/usr/bin/env python -tt
#copyright 2015 Jiong Zhang

"""This file is the main script to call FastLoop 
   fit the loop into the incomplete PDB and write complete PDB
 
   version 0.1 only consider loop in the middle of the structure
"""

print "FastLoop start"

import sys
import os
import shutil

from generateLoopseq import genloopseq

from getTemplate import searchDB
from getTemplate import getTemp

from combineStructure import getLoopfiles
from combineStructure import readParameters

from generateModel import generate_model

FL_HOME = "/Users/jiong/jiong/projects/5_MUFOLD-FL/MUFOLD-FL_py_skMDS"
vmd = "/Applications/LocalApps/VMD1.9.1.app/Contents/MacOS/startup.command"

def main():
   # Get the name from the command line
    if len(sys.argv) == 3 or len(sys.argv) == 4:
    #log = sys.argv[1]
        print 'reading input PDB and sequence file'
        print '......\n'
    else:
        print 'Usage: ${FL_HOME}/FastLoop.py fullseq.fsa target.pdb [savesteps]'
        return

    fullseqfile = sys.argv[1]
    pdbfile = sys.argv[2]
    
    if len(sys.argv) == 4:
        savesteps = sys.argv[3]
    else:
        savesteps = False

    # generate loop sequence file: loopseq.fsa
    genloopseq(fullseqfile, pdbfile)

    # search loop templates
    loopseqfile = "loopseq.fsa"
    loop_length = searchDB(loopseqfile)
    
    fil = open(loopseqfile,'r')
    loopseq = [line for line in fil][-1]
    if not loopseq[-1].isalpha():
        loopseq = loopseq[:-1]
    fil.close()
 
    # frome the output of BLAST copy templates into working directory
    getTemp("blastout.txt", loop_length)
    loopfiles = getLoopfiles()

    # read in full seq file, input pdb file, and loop structure templates
    target, seq, incrd, lpcrds, pdbStart, loopStart, loopEnd = readParameters(fullseqfile, pdbfile, loopfiles)
    print "==============================================================="
    print "Modeling loop for",target,"from residue",loopStart,"to",loopEnd, loopseq
    print "Writing complete structures into folder ./complete_strs"
    print "==============================================================="

    if os.path.exists("./complete_strs"):
        shutil.rmtree("./complete_strs") 
        os.mkdir("./complete_strs")
    else:
        os.mkdir("./complete_strs")

    for loopn in range(len(lpcrds)):
        generate_model(pdbfile, incrd, lpcrds[loopn], loopseq, loopStart, loopEnd, target, loopn, savesteps)
        print "==============================================================="
        print "===              structure optimization done                ==="
        print "==============================================================="
        os.chdir("./complete_strs")
        energymin = vmd+" -dispdev text -e " + FL_HOME + "/FLMD/script/mdrcharm.tcl -args " + target+'complete'+str(loopn)+' '+str(loopStart+pdbStart)+' '+str(loopEnd+pdbStart)
        os.system(energymin)
        os.chdir("../")
    
    # QA 
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
