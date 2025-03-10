#==================================================================================================
# Estimates the per-species correction energy from an evaluation result by performing a linear
# system optimization (least-squares) to minimize excess energy error.
#==================================================================================================

import sys
import os
import time
import argparse
import numpy as np
import ase.io

parser = argparse.ArgumentParser(prog='est-corrections.py')
parser.add_argument('filename', default=None)
args = parser.parse_args()

if args.filename is None:
    print('[ERROR] Filename unspecified; use as "python est-corrections.py <filename>".')
    sys.exit()
elif not os.path.isfile(args.filename):
    print(f'[ERROR] File "{args.filename}" not found!')
    sys.exit()
time_start = time.time()

print(f'Reading "{args.filename}"...')
n_frames = 0
unique_elems = set()
traj = ase.io.iread(args.filename, index=slice(0, None, 1))
for atoms in traj:
    n_frames += 1
    unique_elems.update(set(np.unique(atoms.symbols)))
unique_elems = list(unique_elems)
n_elems = len(unique_elems)
print(f'Found {n_elems} species in {n_frames} frames.')

print(f'Generating count matrix...', end='')
count_mat = np.zeros((n_frames, n_elems), dtype=int)
frame_counter = 0
traj = ase.io.iread(args.filename, index=slice(0, None, 1))
for atoms in traj:
    for i, elem in enumerate(unique_elems):
        count_mat[frame_counter, i] = np.sum(atoms.symbols == elem)
    frame_counter += 1
print('done.')

def is_int(s):
    try:
        int(s)
        return True
    except:
        return False

print(f'Generating error vector...', end='')
error_vec = np.zeros((n_frames,), dtype=float)
frame_counter = 0
with open(args.filename, 'r') as f:
    lines_expected = 0
    header_expected = False
    for line in f:
        if lines_expected == 0 and is_int(line):
            lines_expected = int(line)
            header_expected = True
        elif header_expected:
            header_expected = False
            str_index_start = line.find('REF_energy=') + 11
            str_index_end = line.find(' ', str_index_start)
            true_energy = float(line[str_index_start:str_index_end])
            str_index_start = line.find('MACE_energy=') + 12
            str_index_end = line.find(' ', str_index_start)
            pred_energy = float(line[str_index_start:str_index_end])
            error_vec[frame_counter] = true_energy - pred_energy
            frame_counter += 1
        elif lines_expected > 0:
            pass
print('done.')

print('Calculating corrections:')
corrections_vec = np.linalg.inv(count_mat.T @ count_mat) @ count_mat.T @ error_vec
for i in range(n_elems):
    print(f'    {unique_elems[i]}: {corrections_vec[i]} eV')
print()

time_taken = time.time() - time_start
hours = int(time_taken / 3600)
time_taken -= hours * 3600
mins = int(time_taken / 60)
time_taken -= mins * 60
print(f'Wall time taken (hh:mm:ss) was {hours:02}:{mins:02}:{round(time_taken):02}.')

     