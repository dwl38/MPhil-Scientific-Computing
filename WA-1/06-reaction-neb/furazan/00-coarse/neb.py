#==================================================================================================
# NEB for searching for transition paths at 0 K, e.g. to find activation energy.
#==================================================================================================

import sys
import os
import shutil
import time
import numpy as np
import matplotlib.pyplot as plt
import ase.io

# To suppress warnings from PyTorch coming in through MACE
if not sys.warnoptions:
    import warnings
    warnings.simplefilter('ignore')

from mace.calculators import MACECalculator
from ase.mep.neb import NEB, NEBOptimizer, NEBTools
from ase.optimize import LBFGS
from torch.cuda import is_available as cuda_available
from copy import deepcopy

#==================================================================================================
# Modify parameters here

MAX_IMAGES = 20

ANCHOR_FILES = ['start.xyz', 'anchor_1.xyz', 'anchor_2.xyz', 'anchor_3.xyz',
                'anchor_4.xyz', 'anchor_5.xyz', 'anchor_6.xyz', 'end.xyz']
OPT_DIR = 'opt'
OUTPUT_DIR = 'results'

MACE_MODEL = '../../../00-common/mace-models/MACE-OFF24_medium.model'

#==================================================================================================
# Housekeeping for program launch

time_start = time.time()
DEVICE = 'cuda' if cuda_available() else 'cpu'
print(f'{time.ctime()}: Launched on device "{DEVICE}"')

if os.path.isdir(OUTPUT_DIR):
    shutil.rmtree(OUTPUT_DIR)
    print(f'[WARNING] Pre-existing output directory "{OUTPUT_DIR}" deleted!')
os.mkdir(OUTPUT_DIR)

#==================================================================================================
# Useful functions

def dists_between_images(images):
    lengths = []
    for prev_img, next_img in zip(images, images[1:]):
        lengths.append(np.linalg.norm(next_img.positions - prev_img.positions))
    return np.array(lengths)

def add_intermediate_images(anchors, dist_cutoff, interpolate_method='idpp', max_number=100):
    images = [atoms.copy() for atoms in anchors]
    interp = []
    for i in range(max_number):
        if np.max(dists_between_images(images)) < dist_cutoff:
            break
        dists = dists_between_images(images)
        index = np.argmax(dists)
        to_interpolate = [images[index], images[index].copy(), images[index + 1]]
        neb = NEB(to_interpolate)
        neb.interpolate(method=interpolate_method, apply_constraint=True)
        interp.append([index, to_interpolate[1].copy()])
        images.insert(index + 1, to_interpolate[1].copy())
    return interp, images

#==================================================================================================
# Run NEB

print('Loading MACE-MP-0 model...')
calculator = MACECalculator(model_paths=MACE_MODEL, device=DEVICE, enable_cueq=True)

print('Loading and optimizing anchor states...')
anchors = []
prev_opt_done = os.path.isdir(OPT_DIR)
if prev_opt_done:
    print('  - Found previous results, ', end='')
    for filename in ANCHOR_FILES:
        if not os.path.isfile(os.path.join(OPT_DIR, filename)):
            prev_opt_done = False
    if prev_opt_done:
        print('will load from previously saved results')
        for filename in ANCHOR_FILES:
            atoms = ase.io.read(os.path.join(OPT_DIR, filename))
            atoms.calc = deepcopy(calculator)
            anchors.append(atoms)
    else:
        print('but folder is missing entries - will clear and recalculate')
        shutil.rmtree(OPT_DIR)
if not prev_opt_done:
    os.mkdir(OPT_DIR)
    for i, filename in enumerate(ANCHOR_FILES):
        atoms = ase.io.read(filename)
        atoms.calc = deepcopy(calculator)
        optimizer = LBFGS(atoms, logfile=None)
        optimizer.run(fmax=0.01)
        print(f'  - Optimized {i + 1} out of {len(ANCHOR_FILES)}...')
        ase.io.write(os.path.join(OPT_DIR, filename), atoms)
        anchors.append(atoms)

print('Setting up interpolated states...')
interp, images = add_intermediate_images(anchors, 1e-10, max_number=MAX_IMAGES)
for atoms in images[1:-1]:
    atoms.calc = deepcopy(calculator)
ase.io.write(os.path.join(OUTPUT_DIR, 'interp.xyz'), images)

print('Performing NEB without climb...')
neb = NEB(images=images, climb=False)
optimizer = NEBOptimizer(neb, trajectory=os.path.join(OUTPUT_DIR, 'neb.traj'))
optimizer.run(fmax=1, steps=500)

print('Performing NEB with climb...')
neb.climb = True
optimizer.run(fmax=0.05, steps=700)
ase.io.write(os.path.join(OUTPUT_DIR, 'pathway.xyz'), images)

print('Plotting results...')
fig, axis = plt.subplots()
fig.set_size_inches(12, 9)
nt = NEBTools(images)
nt.plot_band(axis)
fig.savefig(os.path.join(OUTPUT_DIR, 'pathway.png'), dpi=fig.dpi)

nt = NEBTools(ase.io.read(os.path.join(OUTPUT_DIR, 'neb.traj'), ':'))
nt.plot_bands(label=os.path.join(OUTPUT_DIR, 'nebplots'))

#==================================================================================================
# Housekeeping for program end

time_taken = time.time() - time_start
hours = int(time_taken / 3600)
time_taken -= hours * 3600
mins = int(time_taken / 60)
time_taken -= mins * 60
print(f'\n{time.ctime()}: Program ended')
print(f'Wall time taken (hh:mm:ss) was {hours:02}:{mins:02}:{round(time_taken):02}.')
sys.exit()
