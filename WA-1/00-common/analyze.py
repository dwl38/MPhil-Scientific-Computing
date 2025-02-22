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

CONVOLUTION_WIDTH = 100

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

def count_molecules(molecule_dict, target):
    if target[0] == 0:
        if molecule_dict['C'] > 0:
            return 0
        if target[1] == 0:
            if molecule_dict['H'] > 0:
                return 0
            if target[2] == 0:
                if molecule_dict['O'] > 0:
                    return 0
                return int(molecule_dict['N'] / target[3])
            else:
                if molecule_dict['O'] % target[2] != 0:
                    return 0
                possible_count = int(molecule_dict['O'] / target[2])
                return possible_count if molecule_dict['N'] == target[3] * possible_count else 0
        else:
            if molecule_dict['H'] % target[1] != 0:
                return 0
            possible_count = int(molecule_dict['H'] / target[1])
            if molecule_dict['O'] != target[2] * possible_count:
                return 0
            if molecule_dict['N'] != target[3] * possible_count:
                return 0
            return possible_count
    else:
        possible_count = int(molecule_dict['C'] / target[0])
        if molecule_dict['H'] != target[1] * possible_count:
            return 0
        if molecule_dict['O'] != target[2] * possible_count:
            return 0
        if molecule_dict['N'] != target[3] * possible_count:
            return 0
        return possible_count

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
    for target in MOLECULES_TO_FIND:
        composition[target[1]] = 0
    
    molecules = np.unique(molecular_identifiers)
    for molecule_index in molecules:
        molecule_dict = {'C': 0, 'H': 0, 'O': 0, 'N': 0}
        for i in range(N_atoms):
            if molecular_identifiers[i] == molecule_index:
                molecule_dict[atoms[i].symbol] += 1
        for target in MOLECULES_TO_FIND:
            composition[target[1]] += count_molecules(molecule_dict, target[0])
    
    list_of_compositions.append(composition)
    frame_counter += 1

conv_arr = np.linspace(-2, 2, 2 * CONVOLUTION_WIDTH + 1)
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
