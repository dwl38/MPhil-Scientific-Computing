set filename [lindex $argv 0]

# Load file coordinates into VMD
mol new $filename waitfor all

# Open the file to parse lattice vectors
set file [open $filename]
set i 0
while {[gets $file line] >= 0} {
    # Process lines which contain the expression 'Lattice = "...'
    if {[regexp {Lattice=\"[ 0-9\.\-eE]*} $line lattice]} {
        set vectors [regexp -all -inline {\S+} [lindex [split $line \"] 1]]
	set j 0
	foreach v $vectors {
            set v_arr($j) $v
	    incr j
        }
        pbc set "{{$v_arr(0) $v_arr(1) $v_arr(2)} {$v_arr(3) $v_arr(4) $v_arr(5)} {$v_arr(6) $v_arr(7) $v_arr(8)}}" -namd -alignx -first $i -last $i
        incr i
    }
}
close $file
pbc wrap -all

# Graphical settings
set material AOChalky
display projection Orthographic
display depthcue off
axes location Off
display ambientocclusion on
display aoambient 1.000000
display aodirect 0.400000
color change rgb 7 0.000000 0.730000 0.000000
color change rgb 6 0.750000 0.750000 0.750000
color change rgb 2 0.510000 0.510000 0.510000

# Delete default "lines" representation
mol delrep 0 top

# Add CPK representation of atoms
mol representation CPK 1.200000 1.000000 12.000000 12.000000
mol color Name
mol material $material
mol addrep top

