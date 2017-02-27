#calculate rmsd for input str
if {$argc != 2} {
    puts "usage: vmd -dispdev text -e coor2pdb.tcl -args [str]"
    exit
}
	set str [lindex $argv 0]

#    puts $str
    mol new ../common/${str}-allat.psf
    mol addfile ./min-files/minimize-${str}.coor    
    set all [atomselect top all]
    $all writepdb fullatom_${str}.pdb
    exec mv fullatom_${str}.pdb ../
    $all delete
    mol delete top
	
quit


