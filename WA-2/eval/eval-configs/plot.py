#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#==================================================================================================
# Script parameters

INPUT_FILES = ['scratch.xyz',
               'finetune_naive.xyz',
               'finetune_multihead.xyz']

#==================================================================================================
# Do not modify below this line

fig, axes = plt.subplots(2, len(INPUT_FILES))

def is_int(s):
    try:
        int(s)
        return True
    except:
        return False

for i, filename in enumerate(INPUT_FILES):

    ref_energies = list()
    mace_energies = list()
    ref_forces = list()
    mace_forces = list()

    with open(filename, 'r') as f:

        lines_expected = 0
        header_expected = False

        for line in f:

            if lines_expected == 0 and is_int(line):
                lines_expected = int(line)
                header_expected = True
                frame_ref_forces = list()
                frame_mace_forces = list()
            elif header_expected:
                header_expected = False
                str_index_start = line.find('REF_energy=') + 11
                str_index_end = line.find(' ', str_index_start)
                ref_energies.append(float(line[str_index_start:str_index_end]) / lines_expected)
                str_index_start = line.find('MACE_energy=') + 12
                str_index_end = line.find(' ', str_index_start)
                mace_energies.append(float(line[str_index_start:str_index_end]) / lines_expected)
            elif lines_expected > 0:
                lines_expected -= 1
                force_components = line.split()[-6:]
                frame_ref_forces.append([float(force_components[0]), float(force_components[1]), float(force_components[2])])
                frame_mace_forces.append([float(force_components[3]), float(force_components[4]), float(force_components[5])])

            if lines_expected == 0:
                ref_forces.append(np.mean(np.sqrt(np.sum(np.square(frame_ref_forces), axis=-1))))
                mace_forces.append(np.mean(np.sqrt(np.sum(np.square(frame_mace_forces), axis=-1))))

    lower = min(np.min(ref_energies), np.min(mace_energies))
    upper = max(np.max(ref_energies), np.max(mace_energies))
    axes[0][i].plot([lower, upper], [lower, upper], 'k--')
    axes[0][i].plot(ref_energies, mace_energies, 'b.')
    axes[0][i].set_title(f'Energy per atom ({filename.split('.')[0]})')
    axes[0][i].set_xlabel(r'Ref. energy per atom ($eV$)')
    axes[0][i].set_ylabel(r'MACE energy per atom ($eV$)')

    lower = min(np.min(ref_forces), np.min(mace_forces))
    upper = max(np.max(ref_forces), np.max(mace_forces))
    axes[1][i].plot([lower, upper], [lower, upper], 'k--')
    axes[1][i].plot(ref_forces, mace_forces, 'b.')
    axes[1][i].set_title(f'Force ({filename.split('.')[0]})')
    axes[1][i].set_xlabel(r'Ref. force ($eV/\AA$)')
    axes[1][i].set_ylabel(r'MACE force ($eV/\AA$)')

plt.show()

     