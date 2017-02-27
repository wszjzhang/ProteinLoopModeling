
#################################################################
#This is the main script of MDR-amber, run in working directory
#Usage: vmd -dispdev text -e ./mdrcharm.tcl <input pdb>
#################################################################

###############################################################################
#set the path of Rosetta and MDR_package
set maindir /Users/jiong/jiong/projects/5_MUFOLD-FL/MUFOLD-FL_py_skMDS/FLMD/
##############################################################################
if {$argc != 4} {
    puts "usage: vmd -dispdev text -e mdrcharm.tcl -args <pdbin> <loopStart> <loopEnd>"
    exit
}
set pdbin [lindex $argv 0]
set lpStart [lindex $argv 1]
set lpEnd [lindex $argv 2]

proc splitby { string spl_str } {
    set lst [split $string $spl_str]
    for { set cnt 0 } { $cnt < [llength $lst] } { incr cnt } {
        if { [lindex $lst $cnt] == "" } {
            set lst [lreplace $lst $cnt $cnt]
            incr cnt -1
        }
    }
    return $lst
}

# define the procedure that will generate psf and pdb files
proc psf {name dir j loopStart loopEnd } {
    # load the psfgen package
    package require psfgen
    # load the topology file
    topology ../common/top_all27_prot_lipid.inp
    resetpsf

    # change atoms names to match the ones in topology file
    pdbalias atom SER    H HA
    pdbalias atom SER  1HB HB1
    pdbalias atom SER  2HB HB2
    pdbalias atom SER   HG HG1
    pdbalias atom THR    H HN
    pdbalias atom THR 1HG2 HG21
    pdbalias atom THR 2HG2 HG22
    pdbalias atom THR 3HG2 HG23
    pdbalias atom MET    H HN
    pdbalias atom MET  1HB HB1
    pdbalias atom MET  2HB HB2
    pdbalias atom MET  1HG HG1
    pdbalias atom MET  2HG HG2
    pdbalias atom MET  1HE HE1
    pdbalias atom MET  2HE HE2
    pdbalias atom MET  3HE HE3
    pdbalias atom MET   SE SD
    pdbalias atom GLY    H HN
    pdbalias atom GLY  1HA HA1
    pdbalias atom GLY  2HA HA2
    pdbalias atom ALA    H HN
    pdbalias atom ALA  1HB HB1
    pdbalias atom ALA  2HB HB2
    pdbalias atom ALA  3HB HB3
    pdbalias atom ARG    H HN
    pdbalias atom ARG  1HB HB1
    pdbalias atom ARG  2HB HB2
    pdbalias atom ARG  1HG HG1
    pdbalias atom ARG  2HG HG2
    pdbalias atom ARG  1HD HD1
    pdbalias atom ARG  2HD HD2
    pdbalias atom ARG 1HH1 HH11
    pdbalias atom ARG 2HH1 HH12
    pdbalias atom ARG 1HH2 HH21
    pdbalias atom ARG 2HH2 HH22
    pdbalias atom VAL    H HN
    pdbalias atom VAL 1HG1 HG11
    pdbalias atom VAL 2HG1 HG12
    pdbalias atom VAL 3HG1 HG13
    pdbalias atom VAL 1HG2 HG21
    pdbalias atom VAL 2HG2 HG22
    pdbalias atom VAL 3HG2 HG23
    pdbalias atom ILE  CD1 CD
    pdbalias atom ILE    H HN
    pdbalias atom ILE 1HG1 HG11
    pdbalias atom ILE 2HG1 HG12
    pdbalias atom ILE 1HG2 HG21
    pdbalias atom ILE 2HG2 HG22
    pdbalias atom ILE 3HG2 HG23
    pdbalias atom ILE 1HD1 HD1
    pdbalias atom ILE 2HD1 HD2
    pdbalias atom ILE 3HD1 HD3
    pdbalias atom LEU    H HN
    pdbalias atom LEU  1HB HB1
    pdbalias atom LEU  2HB HB2
    pdbalias atom LEU 1HD1 HD11
    pdbalias atom LEU 2HD1 HD12
    pdbalias atom LEU 3HD1 HD13
    pdbalias atom LEU 1HD2 HD21
    pdbalias atom LEU 2HD2 HD22
    pdbalias atom LEU 3HD2 HD23
    pdbalias atom PHE    H HN
    pdbalias atom PHE  1HB HB1
    pdbalias atom PHE  2HB HB2
    pdbalias atom MET   SE S
    
    pdbalias residue HIS HSE
    pdbalias residue MSE MET

    # add the correct N- and C-terminus (if first residue is GLY add a different N-terminus)
    if {[exec head -n 2 ${name} | tail -1 |   cut -b18-20] == "GLY"} {
	segment X {first glyp; last ct3; pdb $name}
    } else {
	segment X {first nter; last ct3; pdb $name}
    }
    coordpdb $name X
    # guess coordinates of hydrogen atoms
    guesscoord
    
    # write the PSF and PDB files
    writepdb ${dir}/${j}-allat.pdb
    writepsf ${dir}/${j}-allat.psf 
    
    # center the pdb file on the origin
    mol new ${dir}/${j}-allat.pdb
    set all [atomselect top all]
    set center [measure center $all]
    set x [expr -1*[lindex $center 0]]
    set y [expr -1*[lindex $center 1]]
    set z [expr -1*[lindex $center 2]]
    
    puts "coor {$x $y $z}"
    $all moveby [vecscale -1 $center]
    $all update
    
    # write final PSF and PDB files 
    set sel [atomselect top all]
    puts [measure center $sel]
    puts [measure minmax $sel]
    $sel writepdb ${dir}/${j}-allat.pdb
    set loop [atomselect top "resid $loopStart to $loopEnd"]    
    $loop set beta 1 
    $sel writepdb ${dir}/fixedrest-${j}.pdb
    $sel delete
    $loop delete
    $all delete 
    mol delete top
}

set dr 0
set names {}; set dirs {};set str {};set cnts {};

# main process start
# create the folders and files needed for the MD simulation
exec pwd
exec mkdir ./common
exec mkdir ${pdbin}-simu

cd ${pdbin}-simu
#exec mkdir heating
exec mkdir simu-files
exec mkdir min-files


#exec cp -pr ${maindir}/MD-select/template/analysis   analysis
exec cp -p  ${maindir}/template/par_all27_prot_lipid.prm ../common
exec cp -p  ${maindir}/template/top_all27_prot_lipid.inp ../common
#set cnts [list 1 2 3 4 5]

# create NAMD configuration files for each simulation and submission scripts for the MD jobs
    
    # create NAMD configuration files for each simulation
    exec echo "set num ${pdbin}" > temp 
    exec cat temp ${maindir}/template/minimize-fixedrest.conf > minimize-fixedrest-${pdbin}.conf

    # create submission scripts for the MD jobs
    exec echo "set name = minimize-fixedrest-${pdbin}" > temp
    exec echo "\$bindir/namd2 +p4 \${name}.conf >& minimize-${pdbin}.out" >> temp

    exec echo "vmd -dispdev text -e ${maindir}/script/coor2pdb.tcl -args ${pdbin}" >> temp

    exec cat ${maindir}/template/job.sh temp > job-${pdbin}.sh
    
    lappend names "../${pdbin}.pdb" 
    lappend dirs  "../common"
    lappend str ${pdbin}


# apply the psf procedure to create all-atom PSF and PDB files
foreach name $names dir $dirs j $str {
    puts "$name $dir $j"
    psf $name $dir $j $lpStart $lpEnd
}

# launch the MD simulations in the queue
#    exec bsub < jobheat-${pdbin}.sh
    exec chmod +x job-${pdbin}.sh
    puts "##############################################################"
    puts "##########       running energy minimization        ##########"
    puts "##############################################################"
    exec ./job-${pdbin}.sh
    cd ..

quit 
