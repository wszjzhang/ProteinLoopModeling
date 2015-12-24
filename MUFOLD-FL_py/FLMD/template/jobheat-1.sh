#!/bin/env tcsh 

#BSUB -a mvapich
#BSUB -J heating
#BSUB -o /dev/null
#BSUB -e /dev/null

#BSUB -n 1

set bindir = /share/apps/bin
