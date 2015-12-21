#!/usr/bin/env python -tt
#copyright 2015 Jiong Zhang

"""
    This script generate fasta file from a pdb file
"""
import sys

def convertAA(s):
    codeAA = {'ALA':'A',
              'ARG':'R',
              'ASN':'N',
              'ASP':'D',
              'ASX':'B',              
              'CYS':'C',              
              'GLU':'E',
              'GLN':'Q',
              'GLX':'Z',
              'GLY':'G',
              'HIS':'H',
              'ILE':'I',
              'LEU':'L',
              'LYS':'K',
              'MET':'M',
              'MSE':'M',
              'PHE':'F',
              'PRO':'P',              
              'SER':'S',
              'THR':'T',
              'TRP':'W',
              'TYR':'Y',
              'PTR':'Y',
              'VAL':'V'}
    return codeAA[s]

def pdb2fasta(pdbfile):
    fil = open(pdbfile,'r')
    ress = [line[17:20] for line in fil if 'CA' in line]
    fil.close()
    
    seq = ""
    for res in ress:
        seq += convertAA(res)

    fil = open("fullseq.fsa", 'w')
    fil.write(">"+pdbfile+"\n")
    fil.write(seq+"\n")
    fil.close()
    
def main():
    pdbfile = sys.argv[1]
    pdb2fasta(pdbfile)

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()


