#!/usr/bin/env python -tt
# Copyright 2015 Jiong Zhang

"""read in the seq file 
   search blastDB
   copy result template in right directory.
"""

import sys
import os
import shutil
import csv

DBPATH = '/Users/jiong/jiong/projects/5_MUFOLD-FL/MUFOLD-FL_py_skMDS/FastLoopDB_NR/'

#input seqfile = 3AA+loop+3AA
def searchDB(seqfile):
    seqf = open(seqfile, 'r')
    for i in range(2):
        seq = seqf.readline()
    loop_length = len(seq)-1
    # search through all databases contain more than aal residues
    print "searching FLDB ..."
    dbs = DBPATH+'blastDB/'+'loop'+str(loop_length)+'AA.fsa'
    for i in range(loop_length, 31):
        dbs = dbs+" "+DBPATH+'blastDB/'+'loop'+str(i)+'AA.fsa'
    blastcmd = 'blastp -db'+'  "'+dbs+'" '+' -query '+seqfile+' -evalue 100 -out blastout.txt'
    os.system(blastcmd)
    print "searching done"
    return loop_length


def getTemp(blastout, match_length):
    #read in blast alginment result
    print "read in alignments from blastout.txt"
    fil = open(blastout, 'r')
    loops = [line.split() for line in fil 
            if '>' in line or 'Query' in line or 'Sbjct' in line][1:]
    fil.close()
    
    #assemble loop candidates with end-to-end distance and start residue number(0 base)
    print "get loop candidates information"
    loop_candidates=[]
    for i in range(0,len(loops),3):
        loopname,loopdis = loops[i][0].split('>')[1].split('|')
        loop_match_query = int(loops[i+1][1])
        loop_match_cand = int(loops[i+2][1])
        loop_start = loop_match_cand - loop_match_query
        if loop_start >= 0:
            loop_candidates.append([loopname,float(loopdis),loop_start])

    for lin in loop_candidates: print lin
    filename = "loop_candidates.txt"
    with open(filename,"wb") as fil:
        writer = csv.writer(fil)
        writer.writerows(loop_candidates)

    #align directory
    try:
        os.stat('templates')
    except:
        os.mkdir('templates')
    #copy loop candidate structure to work directory
    print "copy loop candidate structure to templates directory"
    for loop_cand in loop_candidates:
        loopname = loop_cand[0]
        loop_start = loop_cand[2]
        loop_dis = loop_cand[1]
        folder = loopname[0]
        loop_cand_file = DBPATH+folder+'/'+loopname
        fil = open(loop_cand_file, 'r')
        loop_cand_structure = [line for line in fil][loop_start:loop_start+match_length]
        fil.close()

        #print loopname
        #for lin in loop_cand_structure: print lin
        with open(loopname, 'w') as fil:
            for lin in loop_cand_structure:
                fil.write(lin)
        
        os.rename(loopname,'templates/'+loopname)


def main():
    searchDB(seqfile)    
    getTemp("blastout.txt", loopLength)


if __name__ == '__main__':
    main()


