#==================================================================================================
# Call this script using the following syntax:
#
#     python analyze.py <dirname>
#
# where the dirname can be "./" to run within current directory. This script then attempts to open
# a file named "md.traj" within the working directory (presuming it to be an ASE trajectory file),
# and produces two files in the same directory:
#
#     - traj.xyz, which contains the even-numbered frames of the input trajectory, with coordinates
#       normalized within the periodic boundary conditions;
#
#     - compositions.png, which tracks the molecular compositions over frame number.
#
#==================================================================================================

CUTOFF_LENGTH = 1.8
CUTOFF_LENGTH_SQ = CUTOFF_LENGTH**2

import sys
import os
import numpy as np
import ase.io
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    os.chdir(sys.argv[1])
if not os.path.isfile('md.traj'):
    print('md.traj not found!')
    sys.exit()
if os.path.isfile('traj.xyz'):
    os.remove('traj.xyz')
if os.path.isfile('compositions.png'):
    os.remove('compositions.png')

traj = ase.io.iread('md.traj', index=slice(0, None, 1))

cell_params = None
list_of_compositions = list()
frame_counter = 0

for atoms in traj:

    if cell_params is None:
        cell_params = np.array(atoms.cell.cellpar()[0:3])
    N_atoms = len(atoms)
    
    molecular_identifiers = list(range(N_atoms))
    for i in range(N_atoms):
        for j in range(i + 1, N_atoms):
            dist = atoms.positions[i] - atoms.positions[j]
            dist -= cell_params * np.round(dist / cell_params)
            if np.dot(dist, dist) <= CUTOFF_LENGTH_SQ:
                molecular_identifiers[j] = molecular_identifiers[i]

    composition = dict()
    molecules = np.unique(molecular_identifiers)
    for molecule_index in molecules:
        molecule_contents = {'C': 0, 'H': 0, 'O': 0, 'N': 0}
        for i in range(N_atoms):
            if molecular_identifiers[i] == molecule_index:
                molecule_contents[atoms[i].symbol] += 1
        molecule_name = '$'
        if molecule_contents['C'] > 0:
            molecule_name += 'C_' + str(molecule_contents['C'])
        if molecule_contents['H'] > 0:
            molecule_name += 'H_' + str(molecule_contents['H'])
        if molecule_contents['O'] > 0:
            molecule_name += 'O_' + str(molecule_contents['O'])
        if molecule_contents['N'] > 0:
            molecule_name += 'N_' + str(molecule_contents['N'])
        molecule_name += '$'
        if molecule_name in composition:
            composition[molecule_name] += 1
        else:
            composition[molecule_name] = 1
    list_of_compositions.append(composition)

    if frame_counter % 2 == 0:
        for pos in atoms.positions:
            pos = np.mod(pos, cell_params)
        ase.io.write('traj.xyz', atoms, append=True)
        
    frame_counter += 1

unique_compositions = set()
for composition in list_of_compositions:
    for key in composition.keys():
        unique_compositions.add(key)

fig, axis = plt.subplots()
fig.set_size_inches(12, 9)
for composition in unique_compositions:
    counts = list()
    for frame in list_of_compositions:
        counts.append(frame[composition] if composition in frame else 0)
    axis.plot(counts, label=composition)
axis.set_title('Molecule counts over frame number')
axis.set_xlabel('Frame no.')
axis.set_ylabel('Count')
axis.legend()
fig.savefig('compositions.png', dpi=fig.dpi)

sys.exit()
