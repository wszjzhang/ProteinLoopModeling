#################################################
#    MUFOLD-LOOP                                #
#    organization: University of Missouri       #
#                  MUFOLD team mufold.org       #
#                                               #
#    Version 0.0: Thu Dec 24 15:09:59 CST 2015  #
#    Version 0.1: Sun Sep 18 15:25:54 PDT 2016  #
#    Version 0.1.1: Sun Mar 12 23:11:02 PDT 2017#
#################################################


########################## Reference #################################

Please cite this paper:
MUFOLD-LOOP: a new loop modeling tool for protein structures

########################## Installation ##############################

    This package requires BLAST, numpy, scikit-learn, VMD, NAMD  
    
    Please download the latest version of NAMD from 
        http://www.ks.uiuc.edu/Research/namd/
    
    For BLAST (version 2.2.30+) installation, please refer to 
        https://www.ncbi.nlm.nih.gov/books/NBK279671/

    1, Make sure all the pre-requisite softwares are installed.
    2, In the folder of MUFOLD-LOOP, add the paths of vmd and
       namd in install.py and run ./install.py
    3, Add the path of MUFOLD_LOOP into SHELL enviroment.
    4, Run the command: MUFOLD_LOOP <fullsequence.fasta> <input_pdb>

    Test MUFOLD_LOOP in the test folder, run ./runFL.sh

########################## Usage #####################################

    mufold_loop ./seq.fasta ./str.pdb
    
    You may also refer to the test folder for example of usage. 
    Note: the current version only considers 1 gap in structure.

################# procedure functions description ####################

generateLoopseq.py:
0,  genloopseq(seqfile, pdbfile):
        input file: full seq file, input structure pdb
        generate loop sequence: 3AA+loop+3AA

getTemplate.py:
1,  searchDB(seqfile, loop_length): 
        input seqfile:  3AA+loop+3AA
        search through all databases contain more than aal residues
        output blast result:    blastout.txt

2,  getTemp("blastout.txt", match_length)
        read in blast alignment result
        eliminate incomplete matched sequence
        copy loop candidate structure to work directory

3,  getLoopfiles()
        get loop template files' names

combineStructure.py:
4,  readParameters(seqfile,pdbfile,loopfiles)
        get target,seq,incrd,lpcrds,pdbStart,loopStart,loopEnd        

5,  combStructure(incrd,lpcrds,loopStart,loopEnd)
        for each lpcrd in lpcrds: use procrustes analysis to combine lpcrd into incrd    

6,  convertAA(s):
        convert 1 letter code of AA to 3 letter code
  
7,  writePDB(target, seq, pdbStart, completecrds)
        use seq information and completecrds write pdb

optimizeStructure.py:
8,  pdb2DM(pdbfile)
        read pdb file and calculate the distance matrix

9,  shortest_path(dm)
        reassign adjacent atom pair distances bigger than 3.9 to 3.8
        reassign impropriate pair distances with shortest path 

10, MDSdm2crd(dm)
        use MDS to convert distance matrix to 3D coordinates.


generateModel.py
11, generate_model()
        combine loop with input structure and optimize complete structure
        write out full atom pdb
        call pulchra to generate final structure


FLMD: Energy minimization 
CHARMM force field energy minimization:
    fixed rest part and run energy minimization for 2500 steps


