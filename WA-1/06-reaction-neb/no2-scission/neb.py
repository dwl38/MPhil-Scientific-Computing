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
from ase.neb import NEB, NEBOptimizer, NEBTools
from ase.optimize import LBFGS
from copy import deepcopy

#==================================================================================================
# Modify parameters here

N_IMAGES = 10
START_FILE = 'start.xyz'
END_FILE = 'end.xyz'

OUTPUT_DIR = 'results'

MACE_MODEL = '../00-common/mace-models/MACE-OFF24_medium.model'

#==================================================================================================

time_start = time.time()
if os.path.isdir(OUTPUT_DIR):
    shutil.rmtree(OUTPUT_DIR)
    print(f'[WARNING] Pre-existing output directory "{OUTPUT_DIR}" deleted!')
os.mkdir(OUTPUT_DIR)

print('Loading MACE-MP-0 model...')
calculator = MACECalculator(model_paths=MACE_MODEL, device='cuda', enable_cueq=True)

print('Optimizing start & end states...')
initial = ase.io.read(START_FILE)
initial.calc = deepcopy(calculator)
final = ase.io.read(END_FILE)
final.calc = deepcopy(calculator)
opt_initial = LBFGS(initial)
opt_final = LBFGS(final)
opt_initial.run(fmax=0.01)
opt_final.run(fmax=0.01)

print('Setting up NEB...')
images = [initial,]
images += [initial.copy() for i in range(N_IMAGES - 2)]
images += [final,]
neb = NEB(images=images, climb=False)
neb.interpolate()
for atoms in images[1:-1]:
    atoms.calc = deepcopy(calculator)

print('Performing NEB without climb...')
opt = NEBOptimizer(neb, trajectory=os.path.join(OUTPUT_DIR, 'neb.traj'))
opt.run(fmax=1, steps=500)

print('Performing NEB with climb...')
neb.climb = True
opt.run(fmax=0.05, steps=700)
ase.io.write(os.path.join(OUTPUT_DIR, 'pathway.xyz'), images)

print('Plotting results...')
fig, axis = plt.subplots()
fig.set_size_inches(12, 9)
nt = NEBTools(images)
nt.plot_band(axis)
fig.savefig(os.path.join(OUTPUT_DIR, 'pathway.png'), dpi=fig.dpi)

time_taken = time.time() - time_start
hours = int(time_taken / 3600)
time_taken -= hours * 3600
mins = int(time_taken / 60)
time_taken -= mins * 60
print(f'Wall time taken (hh:mm:ss) was {hours:02}:{mins:02}:{round(time_taken):02}.')