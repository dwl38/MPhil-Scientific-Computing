#!/usr/bin/env python

#==================================================================================================
# Python script for analyzing a trajectory to measure the lattice parameters. Call this script
# using the following syntax:
#
#     python analyze-lat.py <dirname>
#
# where the dirname can be "./" to run within current directory. This script then attempts to open
# a file named "md.traj" within the working directory (presuming it to be an ASE trajectory file),
# and produces two files "lattice_params.png" and "lattice_params.txt" in the same directory, which
# reports the variation of the lattice parameters over time and the final mean values respectively.
#==================================================================================================

import sys
import os
import argparse
import numpy as np
import ase.io
import matplotlib.pyplot as plt

# Parse input arguments
parser = argparse.ArgumentParser(prog='analyze-lat.py')
parser.add_argument('dirname', default='./')
parser.add_argument('--traj', default='md.traj')
parser.add_argument('--index', default=':')
args = parser.parse_args()

# File system handling
print(f'Processing {args.dirname}...', end='')
os.chdir(args.dirname)
if not os.path.isfile(args.traj):
    print(f'{args.traj} not found!')
    sys.exit()
if os.path.isfile('lattice_params.png'):
    os.remove('lattice_params.png')
if os.path.isfile('lattice_params.txt'):
    os.remove('lattice_params.txt')
print(f'found "{args.traj}", writing to "lattice_params.png" and "lattice_params.txt"...')

# Read trajectory
raw_data = list()
traj = ase.io.iread(args.traj, index=args.index)
for atoms in traj:
    raw_data.append(atoms.cell.cellpar())
raw_data = np.array(raw_data)

# Generate plot
fig, axes = plt.subplots(2, 3)
fig.set_size_inches(12, 9)
var_name = ['a', 'b', 'c', r'$\alpha$', r'$\beta$', r'$\gamma$']
y_label_unit = [r'($\AA$)', r'($\degree$)']
for i in range(2):
    for j in range(3):
        axes[i][j].plot(raw_data[:, 3*i+j])
        axes[i][j].set_title('Lattice param ' + var_name[3*i+j] + ' over frame number')
        axes[i][j].set_xlabel('Frame no.')
        axes[i][j].set_ylabel(var_name[3*i+j] + ' ' + y_label_unit[i])
fig.savefig('lattice_params.png', dpi=fig.dpi)

# Generate log
n_data = raw_data.shape[0]
val = [np.mean(raw_data[:, i]) for i in range(6)]
unc = [np.std(raw_data[:, i]) / np.sqrt(n_data - 1) for i in range(6)]
dp = [max(int(np.ceil(-np.log10(u))) + 1, 0) for u in unc]

def sin(a):
    return np.sin(a * np.pi / 180)
def cos(a):
    return np.cos(a * np.pi / 180)

cosines = [cos(angle) for angle in val[3:6]]
cos_product = cosines[0] * cosines[1] * cosines[2]
sines_unc = [(c**2 - cos_product) * sin(a) * (u * np.pi / 180) / c for c, a, u in zip(cosines, val[3:6], unc[3:6])]
prefactor = 1 - np.sum(np.square(cosines)) + (2 * cos_product)
volume = val[0] * val[1] * val[2] * np.sqrt(prefactor)
volume_unc = np.sqrt(((val[1]**2) * (val[2]**2) * prefactor * (unc[0]**2)) +
                     ((val[0]**2) * (val[2]**2) * prefactor * (unc[1]**2)) +
                     ((val[0]**2) * (val[1]**2) * prefactor * (unc[2]**2)) +
                     ((val[0]**2) * (val[1]**2) * (val[2]**2) * (sines_unc[0]**2) / prefactor) +
                     ((val[0]**2) * (val[1]**2) * (val[2]**2) * (sines_unc[1]**2) / prefactor) +
                     ((val[0]**2) * (val[1]**2) * (val[2]**2) * (sines_unc[2]**2) / prefactor))

var_name = ['Length a', 'Length b', 'Length c', 'Angle \u03b1', 'Angle \u03b2', 'Angle \u03b3']
var_unit = [' \u212b', ' \u212b', ' \u212b', '\u00b0', '\u00b0', '\u00b0', ' \u212b\u00b3']
with open('lattice_params.txt', 'w') as output_file:
    for i in range(6):
        output_file.write(f'{var_name[i]} = {round(val[i], dp[i])} \u00b1 {round(unc[i], dp[i])}{var_unit[i]}\n')
    dp = max(int(np.ceil(-np.log10(volume_unc))) + 1, 0)
    output_file.write(f'Volume = {round(volume, dp)} \u00b1 {round(volume_unc, dp)}{var_unit[-1]}\n')

sys.exit()
