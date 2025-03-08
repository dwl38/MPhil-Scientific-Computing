#==================================================================================================
# Measures RDF for water molecules (oxygen atoms)
#==================================================================================================

import sys
import os
import argparse
import numpy as np
import ase.io
import matplotlib.pyplot as plt

# Parse input arguments
parser = argparse.ArgumentParser(prog='analyze.py')
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
if os.path.isfile('rdf.png'):
    os.remove('rdf.png')
if os.path.isfile('rdf.dat'):
    os.remove('rdf.dat')
print(f'found "{args.traj}", writing to "rdf.png" and "rdf.dat"...')

# Read trajectory
dists = list()
traj = ase.io.iread(args.traj, index=args.index)
for atoms in traj:
    lattice_params = np.array(atoms.cell.cellpar()[0:3], dtype=float)
    oxygens = atoms.positions[atoms.symbols == 'O']
    for pos in oxygens:
        disp = oxygens - pos
        disp -= lattice_params * np.round(disp / lattice_params)
        dists += list(np.linalg.norm(disp, axis=-1))
dists = np.array(dists)[np.nonzero(dists)]

# Calculate RDFs
max_dist = np.max(dists)
n_dists = dists.shape[0]
bin_edges = np.linspace(0.0, max_dist, 101)
left_bin_edges = bin_edges[:-1]
right_bin_edges = bin_edges[1:]
bin_centers = (left_bin_edges + right_bin_edges) / 2
rdf = np.zeros((100,), dtype=float)
for i in range(100):
    count = np.sum(np.logical_and(dists > left_bin_edges[i], dists <= right_bin_edges[i]))
    rdf[i] = 3 * count / (4 * np.pi * n_dists * ((right_bin_edges[i]**3) - (left_bin_edges[i]**3)))

# Generate log
with open('rdf.dat', 'w') as output_file:
    output_file.write(f'# RDF for {args.dirname}\n')
    output_file.write('# r g(r)\n')
    for r, g in zip(bin_centers, rdf):
        output_file.write(f'{r} {g}\n')

# Generate plot
fig, axis = plt.subplots()
fig.set_size_inches(12, 9)
axis.plot(np.concatenate(([0,], bin_centers)), np.concatenate(([0,], rdf)))
axis.set_title(f'RDF for {args.dirname}')
axis.set_xlabel(r'r ($\AA$)')
axis.set_ylabel(r'$g(r)$ ($\AA^{-3}$)')
fig.savefig('rdf.png', dpi=fig.dpi)

sys.exit()
