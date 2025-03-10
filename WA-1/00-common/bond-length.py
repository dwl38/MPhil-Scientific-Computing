#!/usr/bin/env python

#==================================================================================================
# Short Python script which analyzes a trajectory, and plots the distribution of bond lengths
# (flattened over time) for all pairs of elements.
#==================================================================================================

import sys
import os
import time
import argparse
import numpy as np
import ase.io
import matplotlib.pyplot as plt

# Parse input arguments
parser = argparse.ArgumentParser(prog='bond-length.py')
parser.add_argument('dirname', default='./')
parser.add_argument('--traj', default='md.traj')
parser.add_argument('-o', '--output', default='bond_lengths.png')
args = parser.parse_args()

# File system handling
os.chdir(args.dirname)
print(f'Running in {args.dirname}...', end='')
if not os.path.isfile(args.traj):
    print(f'{args.traj} not found!')
    sys.exit()
if os.path.isfile(args.output):
    os.remove(args.output)

# Read trajectory
time_start = time.time()
lengths = {'C': {'C': [], 'H': [], 'O': [], 'N': []},
           'H': {'C': [], 'H': [], 'O': [], 'N': []},
           'O': {'C': [], 'H': [], 'O': [], 'N': []},
           'N': {'C': [], 'H': [], 'O': [], 'N': []}}
traj = ase.io.iread(args.traj, index=slice(0, None, 1))
for atoms in traj:
    cell_params = np.array(atoms.cell.cellpar()[0:3])
    N_atoms = len(atoms)
    for i in range(N_atoms):
        dist = atoms.positions - atoms.positions[i]
        dist -= cell_params * np.round(dist / cell_params)
        dist = np.sqrt(np.sum(np.square(dist), axis=1))
        for j in range(i + 1, N_atoms):
            lengths[atoms[i].symbol][atoms[j].symbol].append(dist[j])

# Symmetricize lengths matrix
elems = ['C', 'H', 'O', 'N']
for i in range(4):
    for j in range(i + 1, 4):
        lengths[elems[i]][elems[j]] += lengths[elems[j]][elems[i]]
        lengths[elems[j]][elems[i]] = None

# Plotting function
def plot_histogram(axis, lengths, elem_i, elem_j):
    samples = np.array(lengths[elem_i][elem_j])
    if samples.shape[0] > 3300000:
        np.random.shuffle(samples)
        samples = samples[0:3300000]
    lower = np.min(samples)
    upper = np.max(samples)
    width = (upper - lower) / np.sqrt(samples.shape[0])
    scalefactor = -0.5 * np.power(width, -2)
    prefactor = 1.0 / (np.sqrt(2 * np.pi) * width * samples.shape[0])
    x = np.linspace(lower, upper, 250)
    y = prefactor * np.sum(np.exp(scalefactor * np.square(x[:, np.newaxis] - samples[np.newaxis, :])), axis=-1)
    height = np.max(y)
    peak = x[np.argmax(y)]
    axis.plot(x, y, 'b-', label=f'{elem_i}-{elem_j}')
    axis.plot([peak, peak], [0, height], '--', color='gray', label=f'{round(peak, 2)} ' + r'$\AA$')
    axis.set_title(f'Bond length distributions for {elem_i}-{elem_j}')
    axis.set_xlabel(r'Bond length ($\AA$)')
    axis.set_ylabel('Density')
    axis.legend()

# Generate plot
fig, axes = plt.subplots(10, 1)
fig.set_size_inches(12, 70)
counter = 0
for i in range(4):
    for j in range(i, 4):
        plot_histogram(axes[counter], lengths, elems[i], elems[j])
        counter += 1
fig.savefig(args.output, dpi=fig.dpi)

# End of program
time_taken = time.time() - time_start
hours = int(time_taken / 3600)
time_taken -= hours * 3600
mins = int(time_taken / 60)
time_taken -= mins * 60
print(f'wall time taken (hh:mm:ss) was {hours:02}:{mins:02}:{round(time_taken):02}.')

sys.exit()

