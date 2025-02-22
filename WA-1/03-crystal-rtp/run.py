#==================================================================================================
# MD simulation at 300 K, 1 atm over 100 ps to demonstrate stability of solid crystal.
#==================================================================================================

import sys
import time
import ase.io

# To suppress warnings from PyTorch coming in through MACE
if not sys.warnoptions:
    import warnings
    warnings.simplefilter('ignore')

from mace.calculators import MACECalculator
from ase.md import MDLogger
from ase.md.npt import NPT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units

#==================================================================================================
# Modify parameters here

N_STEPS = 100000
TIMESTEP = 1 * units.fs
START_FILE = 'start.xyz'

MACE_MODEL = '../00-common/mace-models/MACE-OFF24_medium.model'

#==================================================================================================
# Script

print('Loading MACE-MP-0 model...')
calculator = MACECalculator(model_paths=MACE_MODEL, device='cuda', enable_cueq=True)
atoms = ase.io.read(START_FILE)
atoms.calc = calculator

print('Starting MD...')
time_start = time.time()
MaxwellBoltzmannDistribution(atoms, temperature_K=300)
dyn = NPT(atoms, timestep=TIMESTEP, temperature_K=300, externalstress=1.01325 * units.bar,
          ttime=25 * units.fs, pfactor=((75 * units.fs)**2) * (100 * units.GPa))
dyn.attach(MDLogger(dyn, atoms, 'md.log', header=True, stress=True, peratom=True, mode='w'),
           interval=100)
dyn.attach(ase.io.Trajectory('md.traj', 'w', atoms).write, interval=10)

print('Running MD...')
for i in range(20):
    print(f'    Progress: {int(i * 5)}%...')
    dyn.run(int(N_STEPS/20))
    ase.io.write('final.xyz', atoms, append=False)
print('MD finished!')

time_taken = time.time() - time_start
hours = int(time_taken / 3600)
time_taken -= hours * 3600
mins = int(time_taken / 60)
time_taken -= mins * 60
print(f'Wall time taken (hh:mm:ss) was {hours:02}:{mins:02}:{round(time_taken):02}.')

