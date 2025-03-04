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
# Internally, the script represents molecules using a graph-theoretic representation, with
# abbreviated labels for common substructures.
#
#==================================================================================================

import sys
import os
import time
import argparse
import numpy as np
import ase.io
import matplotlib.pyplot as plt


#==================================================================================================
# Script parameters

CUTOFF_LENGTHS = {'C': {'C': 1.80, 'H': 1.36, 'O': 1.79, 'N': 1.90},
                  'H': {'C': 1.36, 'H': 0.93, 'O': 1.20, 'N': 1.26},
                  'O': {'C': 1.79, 'H': 1.20, 'O': 1.51, 'N': 1.56},
                  'N': {'C': 1.90, 'H': 1.26, 'O': 1.56, 'N': 1.38}}

CONVOLUTION_WIDTH = 300

MOLECULES_TO_FIND = [[[6, 6, 6, 6], r'TATB', 'TATB'],       # TATB, C6H6O6N6
                     [[6, 4, 5, 6], r'TATB-F1', 'TATB-F1'], # Furazan intermediate, C6H4O5N6
                     [[6, 2, 4, 6], r'TATB-F2', 'TATB-F2'], # Furazan intermediate, C6H2O4N6
                     [[6, 0, 3, 6], r'TATB-F3', 'TATB-F3'], # Furazan intermediate, C6O3N6
                     [[1, 0, 2, 0], r'$CO_{2}$', 'CO2'],    # Carbon dioxide
                     [[1, 0, 1, 0], r'$CO$', 'CO'],         # Carbon monoxide
                     [[0, 2, 1, 0], r'$H_{2}O$', 'H2O'],    # Water
                     [[0, 0, 2, 1], r'$NO_{2}$', 'NO2'],    # Nitrous dioxide
                     [[0, 0, 1, 2], r'$N_{2}O$', 'N2O'],    # Nitrous oxide
                     [[0, 0, 1, 1], r'$NO$', 'NO'],         # Nitrogen monoxide
                     [[0, 0, 0, 2], r'$N_{2}$', 'N2']]      # Nitrogen gas


#==================================================================================================

# Parse input arguments
parser = argparse.ArgumentParser(prog='analyze.py')
parser.add_argument('dirname', default='./')
parser.add_argument('--traj', default='md.traj')
parser.add_argument('-o', '--output', default='compositions.png')
parser.add_argument('--filter', action='extend', nargs='*', default=None)
parser.add_argument('--conv-width', type=int, default=CONVOLUTION_WIDTH)
args = parser.parse_args()

# File system handling
print(f'Processing {args.dirname}...', end='')
os.chdir(args.dirname)
if not os.path.isfile(args.traj):
    print(f'{args.traj} not found!')
    sys.exit()
if os.path.isfile(args.output):
    os.remove(args.output)
print(f'found {args.traj}, writing to {args.output}...')

# Filter only desired molecules
if args.filter is not None:
    whitelist = [[np.array(mol[0]), np.nonzero(mol[0])[0], mol[1]] for mol in MOLECULES_TO_FIND if mol[2] in args.filter]
else:
    whitelist = [[np.array(mol[0]), np.nonzero(mol[0])[0], mol[1]] for mol in MOLECULES_TO_FIND]
MOLECULES_TO_FIND = whitelist
print('Filtering for only: ', end='')
for mol in MOLECULES_TO_FIND:
    print(mol[2], end=' ')
print()

# Internal representation of adjacency matrix
elems = CUTOFF_LENGTHS.keys()
CUTOFF_LENGTHS_SQ = np.zeros((4, 4), dtype=float)
for i in range(4):
    for j in rang(4):
        CUTOFF_LENGTHS_SQ[i, j] = (CUTOFF_LENGTHS[elems[i]][elems[j]])**2
CONVOLUTION_WIDTH = args.conv_width

# Read trajectory
time_start = time.time()
list_of_compositions = list()
frame_counter = 0
molecule = np.zeros((4,), dtype=int)
traj = ase.io.iread(args.traj, index=slice(0, None, 1))

for atoms in traj:

    # Get parameters for this frame
    cell_params = np.array(atoms.cell.cellpar()[0:3])
    N_atoms = len(atoms)
    
    # Internal representation
    atoms.numbers[atoms.symbols == 'C'] = 0
    atoms.numbers[atoms.symbols == 'H'] = 1
    atoms.numbers[atoms.symbols == 'O'] = 2
    atoms.numbers[atoms.symbols == 'N'] = 3
    
    # Find connected components
    molecular_label = np.linspace(0, N_atoms, N_atoms, dtype=int)
    for i in range(N_atoms):
        dist = atoms.positions - atoms.positions[i]
        dist -= cell_params * np.round(dist / cell_params)
        dist_sq = np.sum(np.square(dist), axis=1)
        molecular_label[dist_sq < CUTOFF_LENGTHS_SQ[atoms.numbers[i]]] = molecular_label[i]

    # Deduce compositions based on largest fingerprints
    composition = {key[2]: 0 for key in MOLECULES_TO_FIND}
    unique_labels = np.unique(molecular_label)
    for index in unique_labels:
        molecule[:] = 0
        molecule[0] += np.sum(np.logical_and(molecular_label == index, atoms.numbers == 0))
        molecule[1] += np.sum(np.logical_and(molecular_label == index, atoms.numbers == 1))
        molecule[2] += np.sum(np.logical_and(molecular_label == index, atoms.numbers == 2))
        molecule[3] += np.sum(np.logical_and(molecular_label == index, atoms.numbers == 3))
        for target in MOLECULES_TO_FIND:
            molecule_components = molecule[target[1]]
            count = np.min(np.floor_divide(molecule_components, target[0][target[1]]))
            composition[target[2]] += count
            molecule -= count * target[0]
    
    # Add to found compositions and continue iterating
    list_of_compositions.append(composition)
    frame_counter += 1

# Convolution array for Gaussian smoothing
conv_arr = np.linspace(-3, 3, 3 * CONVOLUTION_WIDTH + 1)
conv_arr = np.exp(-0.5 * np.square(conv_arr))
conv_arr = conv_arr / np.sum(conv_arr)

# Generate plot
fig, axis = plt.subplots()
fig.set_size_inches(12, 9)
for target in MOLECULES_TO_FIND:
    counts = [frame[target[2]] for frame in list_of_compositions]
    if np.sum(counts) > 0:
        axis.plot(np.convolve(counts, conv_arr, 'valid'), label=target[2])
axis.set_title('Molecule counts over frame number')
axis.set_xlabel('Frame no.')
axis.set_ylabel('Count')
axis.legend()
fig.savefig(args.output, dpi=fig.dpi)

# End of program
time_taken = time.time() - time_start
hours = int(time_taken / 3600)
time_taken -= hours * 3600
mins = int(time_taken / 60)
time_taken -= mins * 60
print(f'Wall time taken (hh:mm:ss) was {hours:02}:{mins:02}:{round(time_taken):02}.')

sys.exit()
