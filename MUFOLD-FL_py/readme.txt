procedure functions description

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
        reassigne adjacent atom pair distances bigger than 3.9 to 3.8
        reassign impropriate pair distances with shortest path 

10, MDSdm2crd(dm)
        use MDS to convert distance matrix to 3D coordinates.


generateModel.py
11, generate_model()
        combine loop with input structure and optimize complete structure
        write out full atom pdb
        call pulchra to generate final structure
