#!/usr/bin/env python -tt
from __future__ import print_function

import sys, stat
import os
from tempfile import mkstemp
from shutil import move
from os import remove, close

from sys import platform

############################## configuration #########################

vmd_path = "/Applications/LocalApps/VMD\ 1.9.app/Contents/MacOS/startup.command"
namd_path = "/Users/jiongz/jiong/projects/5_MUFOLD-FL/NAMD_2.12_MacOSX-x86_64-multicore"

######################################################################

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
    if len(sys.argv) != 1:
        print("install: ./install.py ")
        return

    if platform == "linux" or platform == "linux2":
        # linux
        os_version = 'linux'
    elif platform == "darwin":
        # OS X
        os_version = 'osx'
    elif platform == "win32":
        # Windows...
        os_version = 'win32'


    cwd = os.getcwd()
    replace("./structure_generation/generateModel.py","PULCHRA = '", "    PULCHRA = '"+cwd+"/pulchra304/bin/"+os_version+"/pulchra '\n")
    replace("./structure_generation/getTemplate.py","DBPATH = ", "DBPATH = '"+cwd+"/FastLoopDB_NR/'\n")

    replace("./mufold_loop","FL_HOME = ", "FL_HOME = '"+cwd+"'\n")
    replace("./mufold_loop","vmd = ", "vmd = '"+vmd_path+"'\n")

    replace("./FLMD/script/mdrundir.sh", "set FL_HOME", "set FL_HOME = '"+cwd+"'\n")
    replace("./FLMD/script/mdrcharm.tcl", "set maindir", 'set maindir "'+cwd+'/FLMD\n')
    replace("./FLMD/template/job.sh", "set bindir", 'set bindir = "'+namd_path+'\n')

    os.chmod("./mufold_loop", stat.S_IXUSR | stat.S_IRUSR | stat.S_IWUSR)
    print("Wrote paths in scripts")
    print("Please make sure blast is installed")

if __name__ == '__main__':
    main()
