#!/usr/bin/env python

#==================================================================================================
# Python script for analyzing a trajectory, and identifiying the molecular species present in each
# frame. Call this script using the following syntax:
#
#     python analyze.py <dirname>
#
# where the dirname can be "./" to run within current directory. This script then attempts to open
# a file named "md.traj" within the working directory (presuming it to be an ASE trajectory file),
# and produces a file "compositions.png" in the same directory, which tracks the molecular
# compositions over frame number.
#
#--------------------------------------------------------------------------------------------------
#
# Internally, the script represents molecules using a list of four integers [C, H, O, N],
# representing the number of carbon, hydrogen, oxygen, and nitrogen atoms in the molecule.
#
#==================================================================================================

import sys
import os
import numpy as np
import ase.io
import matplotlib.pyplot as plt


#==================================================================================================
# Script parameters

CUTOFF_LENGTH = 1.85
CUTOFF_LENGTH_SQ = CUTOFF_LENGTH**2

CONVOLUTION_WIDTH = 200

MOLECULES_TO_FIND = [[[6, 6, 6, 6], r'TATB'],     # TATB, C6H6O6N6
                     [[6, 4, 5, 6], r'TATB-F1'],  # TATB-F1 (furazan intermediate), C6H4O5N6
                     [[6, 2, 4, 6], r'TATB-F2'],  # TATB-F2 (furazan intermediate), C6H2O4N6
                     [[6, 0, 3, 6], r'TATB-F3'],  # TATB-F3 (furazan intermediate), C6O3N6
                     [[1, 0, 2, 0], r'$CO_{2}$'], # Carbon dioxide, CO2
                     [[1, 0, 1, 0], r'$CO$'],     # Carbon monoxide, CO
                     [[0, 2, 1, 0], r'$H_{2}O$'], # Water, H2O
                     [[0, 0, 0, 2], r'$N_{2}$'],  # Nitrogen gas, N2
                     [[0, 0, 1, 2], r'$N_{2}O$'], # Nitrous oxide, N2O
                     [[0, 0, 1, 1], r'$NO$'],     # Nitrogen monoxide, NO
                     [[0, 0, 2, 1], r'$NO_{2}$']] # Nitrogen dioxide, NO2


#==================================================================================================

if len(sys.argv) > 1:
    os.chdir(sys.argv[1])
if not os.path.isfile('md.traj'):
    print('md.traj not found!')
    sys.exit()
if os.path.isfile('compositions.png'):
    os.remove('compositions.png')

def count_molecules(molecule, target):
    if np.sum(np.mod(molecule, target)) > 0:
        return 0
    nonzero_index = np.nonzero(target)[0][0]
    possible_count = molecule[nonzero_index] / target[nonzero_index]
    for i in range(4):
        if target[i] * possible_count != molecule[i]:
            return 0
    return possible_count

traj = ase.io.iread('md.traj', index=slice(0, None, 1))

list_of_compositions = list()
frame_counter = 0

for atoms in traj:

    cell_params = np.array(atoms.cell.cellpar()[0:3])
    N_atoms = len(atoms)
    
    molecular_identifiers = np.linspace(0, N_atoms, N_atoms, dtype=int)
    for i in range(N_atoms):
        dist = atoms.positions - atoms.positions[i]
        dist -= cell_params * np.round(dist / cell_params)
        dist_sq = np.sum(np.square(dist), axis=1)
        molecular_identifiers[dist_sq <= CUTOFF_LENGTH_SQ] = molecular_identifiers[i]

    composition = {key[1]: 0 for key in MOLECULES_TO_FIND}
    
    molecules = np.unique(molecular_identifiers)
    for molecule_index in molecules:
        molecule_contents = [0, 0, 0, 0]
        molecule_contents[0] += np.sum(np.logical_and(molecular_identifiers == molecules, atoms.symbols == 'C'))
        molecule_contents[1] += np.sum(np.logical_and(molecular_identifiers == molecules, atoms.symbols == 'H'))
        molecule_contents[2] += np.sum(np.logical_and(molecular_identifiers == molecules, atoms.symbols == 'O'))
        molecule_contents[3] += np.sum(np.logical_and(molecular_identifiers == molecules, atoms.symbols == 'N'))
        for target in MOLECULES_TO_FIND:
            composition[target[1]] += count_molecules(molecule_contents, target[0])
    
    list_of_compositions.append(composition)
    frame_counter += 1

conv_arr = np.linspace(-3, 3, 3 * CONVOLUTION_WIDTH + 1)
conv_arr = np.exp(-0.5 * np.square(conv_arr))
conv_arr = conv_arr / np.sum(conv_arr)

fig, axis = plt.subplots()
fig.set_size_inches(12, 9)
for target in MOLECULES_TO_FIND:
    counts = list()
    display = False
    for frame in list_of_compositions:
        counts.append(frame[target[1]])
        if frame[target[1]] > 0:
            display = True
    if display:
        axis.plot(np.convolve(counts, conv_arr, 'valid'), label=target[1])
axis.set_title('Molecule counts over frame number')
axis.set_xlabel('Frame no.')
axis.set_ylabel('Count')
axis.legend()
fig.savefig('compositions.png', dpi=fig.dpi)

sys.exit()
