# Use with VMD in the format:
#   vmd <filename> -e view.tcl

set material AOChalky
display backgroundgradient on
color Display Background white
display projection Orthographic
display depthcue off
axes location Off
display ambientocclusion on
display aoambient 1.000000
display aodirect 0.400000
color change rgb 7 0.000000 0.730000 0.000000
color change rgb 6 0.750000 0.750000 0.750000
color change rgb 2 0.510000 0.510000 0.510000

# Ensure wrapping at periodic boundary conditions
pbc set {13.9116 13.9116 13.9116} -all
pbc wrap all

# Delete default "lines" representation
mol delrep 0 top

# Change color of hydrogen & carbon atoms
color Name H silver
color Name C black

# Add CPK representation of atoms
mol representation CPK 1.20000 0.000000 12.000000 12.000000
mol color Name
mol material $material
mol addrep top

# Add dynamic representation of bonds
mol representation DynamicBonds 1.700000 0.100000 12.000000
mol color Name
mol material $material
mol addrep top

