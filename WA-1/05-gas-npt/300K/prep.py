import sys
import time
import numpy as np
import ase.io

# To suppress warnings from PyTorch coming in through MACE
if not sys.warnoptions:
    import warnings
    warnings.simplefilter('ignore')

from mace.calculators import MACECalculator
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import MDLogger
from ase.md.npt import NPT
from ase import units

#==================================================================================================
# Modify parameters here

N_STEPS = 5000
TIMESTEP = 1 * units.fs

MACE_MODEL = '../../00-common/mace-models/MACE-OFF24_medium.model'
START_FILE = '../../00-common/geoms/crystal.xyz'

#==================================================================================================
# Script

print('Preparing system...')
atoms = ase.io.read(START_FILE)
width_x = np.max(atoms.positions[:,0]) - np.min(atoms.positions[:,0]) + 1
width_y = np.max(atoms.positions[:,1]) - np.min(atoms.positions[:,1]) + 1
width_z = np.max(atoms.positions[:,2]) - np.min(atoms.positions[:,2]) + 1
width = max(max(width_x, width_y), width_z)
atoms.set_cell(width * np.identity(3))
atoms.set_pbc(True)

print('Loading MACE-MP-0 model...')
calculator = MACECalculator(model_paths=MACE_MODEL, device='cuda', enable_cueq=True)
atoms.calc = calculator

print('Starting MD...')
time_start = time.time()
MaxwellBoltzmannDistribution(atoms, temperature_K=300)
dyn = NPT(atoms, timestep=TIMESTEP, temperature_K=300, externalstress=1.01325 * units.bar,
          ttime=25 * units.fs, pfactor=((75 * units.fs)**2) * (100 * units.GPa))

print('Running MD...')
dyn.run(N_STEPS)
ase.io.write('final.xyz', atoms, append=False)
print('MD finished!')

time_taken = time.time() - time_start
hours = int(time_taken / 3600)
time_taken -= hours * 3600
mins = int(time_taken / 60)
time_taken -= mins * 60
print(f'Wall time taken (hh:mm:ss) was {hours:02}:{mins:02}:{round(time_taken):02}.')
