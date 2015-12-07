#!/usr/bin/env python -tt

import os
from tempfile import mkstemp
from shutil import move
from os import remove, close

def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with open(abs_path,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                if pattern in line:
                    new_file.write(line.replace(line, subst))
                else:
                    new_file.write(line)
    close(fh)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

def main():
    cwd = os.getcwd()
    replace("./generateModel.py","PULCHRA = '", "    PULCHRA = '"+cwd+"/pulchra304/bin/osx/pulchra '\n")
    replace("./getTemplate.py","DBPATH = ", "DBPATH = '"+cwd+"/FastLoopDB_NR/'\n")
    print "Wrote paths in scripts"
    print "Please make sure blast is installed"

if __name__ == '__main__':
    main()
