#calculate rmsd for input str
if {$argc != 2} {
    puts "usage: vmd -dispdev text -e single-avrmsd.tcl -args <pdbin> "
    exit
}
	set str [lindex $argv 0]

#    puts $str
    mol new ../common/${str}-allat.pdb
    animate read dcd ../${str}-simu/heating-files/heating-${str}.dcd
    set m1 [atomselect top "name CA" frame 0]
    set m2 [atomselect top "name CA"]
    set nf [molinfo top get numframes]
    set avrmsd 0
    for {set i 1} {$i <= $nf} {incr i} {
        $m2 frame $i; $m2 update
        $m2 move [measure fit $m2 $m1]
        set rmsd [measure rmsd $m1 $m2]
        set avrmsd [expr $avrmsd+$rmsd]
    }
    puts $nf
    set avrmsd [expr $avrmsd/$nf]
    $m1 delete
    $m2 delete
    mol delete top
	
	exec echo $str $avrmsd >> ../c-avrmsd.dat
quit


