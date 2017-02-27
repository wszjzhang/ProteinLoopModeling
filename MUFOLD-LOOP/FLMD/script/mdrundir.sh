#!/bin/tcsh

#run in work dir

set strs = `ls -1 *.pdb | cut -d"." -f1`
set FL_HOME = '/Users/jiong/jiong/projects/5_MUFOLD-FL/MUFOLD-FL_py_skMDS'

foreach str ( $strs )
   vmd -dispdev text -e ${FL_HOME}/FLMD/script/mdrcharm.tcl -args $str $1 $2
end
