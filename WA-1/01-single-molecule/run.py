#==================================================================================================
# MD simulation at user-specified temperature (NVT ensemble) of one molecule over 1 ns (1000 ps)
#==================================================================================================

import sys
import os
import time
import ase.io

# To suppress warnings from PyTorch coming in through MACE
if not sys.warnoptions:
    import warnings
    warnings.simplefilter('ignore')

from mace.calculators import MACECalculator
from ase.md import MDLogger
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units

#==================================================================================================
# Modify parameters here

N_STEPS = 1000000
TIMESTEP = 1 * units.fs
START_FILE = 'start.xyz'

MACE_MODEL = '../00-common/mace-models/MACE-OFF24_medium.model'

#==================================================================================================
# Script

print('Loading MACE-MP-0 model...')
calculator = MACECalculator(model_paths=MACE_MODEL, device='cuda', enable_cueq=True)
atoms = ase.io.read(START_FILE)
atoms.calc = calculator
atoms.set_pbc(True)

if len(sys.argv) > 1:
    temp = float(sys.argv[1])
else:
    temp = input('\nInput temperature in K: ')
    temp = float(temp)

if temp < 0:
    raise RuntimeError('Invalid temperature!')
folder_prefix = str(int(temp)) + 'K'

if(os.path.isdir(folder_prefix)):
    folder_prefix += '_1'
    attempt = 2
    while(os.path.isdir(folder_prefix)):
        folder_prefix[:-1] += str(attempt)
        attempt += 1
os.mkdir(folder_prefix)

print('Starting MD...')
time_start = time.time()
MaxwellBoltzmannDistribution(atoms, temperature_K=temp)
dyn = Langevin(atoms, timestep=TIMESTEP, temperature_K=temp, friction=(0.01/TIMESTEP))
dyn.attach(MDLogger(dyn, atoms, os.path.join(folder_prefix, 'md.log'),
                    header=True, stress=False, peratom=True, mode='w'), interval=100)
traj = ase.io.Trajectory(os.path.join(folder_prefix, 'md.traj'), 'w', atoms)
dyn.attach(traj.write, interval=10)

print('Running MD...')
for i in range(20):
    print(f'    Progress: {int(i * 5)}%...')
    dyn.run(int(N_STEPS/20))
print('MD finished!')
ase.io.write(os.path.join(folder_prefix, 'final.xyz'), atoms)

time_taken = time.time() - time_start
hours = int(time_taken / 3600)
time_taken -= hours * 3600
mins = int(time_taken / 60)
time_taken -= mins * 60
print(f'Wall time taken (hh:mm:ss) was {hours:02}:{mins:02}:{round(time_taken):02}.')

