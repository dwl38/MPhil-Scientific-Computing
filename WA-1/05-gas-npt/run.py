#==================================================================================================
# MD simulation at 300 K, 1 atm over 100 ps to demonstrate stability of solid crystal.
#==================================================================================================

import sys
import os
import shutil
import time
import ase.io

# To suppress warnings from PyTorch coming in through MACE
if not sys.warnoptions:
    import warnings
    warnings.simplefilter('ignore')

from mace.calculators import MACECalculator
from ase.md import MDLogger
from ase.md.npt import NPT
from ase import units

#==================================================================================================
# Modify parameters here

TIMESTEP = 1 * units.fs

MACE_MODEL = '../00-common/mace-models/MACE-OFF24_medium.model'

#==================================================================================================
# Script

start_temp = int(sys.argv[1]) if len(sys.argv) > 1 else 300
end_temp = int(sys.argv[2]) if len(sys.argv) > 2 else 400
steps_per_degree = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
equilibration_steps = int(sys.argv[4]) if len(sys.argv) > 4 else 0
print(f'Running crystal gradual heat from {start_temp}K to {end_temp}K with {steps_per_degree}' +
      f' timesteps per degree and {equilibration_steps} final equilibration timesteps...')

start_folder = str(start_temp) + 'K'
end_folder = str(end_temp) + 'K'
if not os.path.isdir(start_folder) or not os.path.isfile(os.path.join(start_folder, 'final.xyz')):
    print(f'[ERROR] Start folder {start_folder} does not exist!')
    raise RuntimeError(f'[ERROR] Start folder "{start_folder}" does not exist!')
    sys.exit()
if os.path.isdir(end_folder):
    shutil.rmtree(end_folder)
    print(f'[WARNING] Pre-existing end folder "{end_folder}" deleted!')
os.mkdir(end_folder)

print('Loading MACE-MP-0 model...')
calculator = MACECalculator(model_paths=MACE_MODEL, device='cuda', enable_cueq=True)
atoms = ase.io.read(os.path.join(start_folder, 'final.xyz'))
atoms.calc = calculator

print('Starting MD...')
time_start = time.time()
dyn = NPT(atoms, timestep=TIMESTEP, temperature_K=start_temp, externalstress=1.01325 * units.bar,
          ttime=25 * units.fs, pfactor=((75 * units.fs)**2) * (100 * units.GPa))
dyn.attach(MDLogger(dyn, atoms, os.path.join(end_folder, 'md.log'),
           header=True, stress=True, peratom=True, mode='w'), interval=100)
dyn.attach(ase.io.Trajectory(os.path.join(end_folder, 'md.traj'), 'w', atoms).write, interval=10)

print('Running MD...')
for t in range(start_temp, end_temp, 1):
    print(f'    Progress: running {t+1}K...')
    dyn.set_temperature(temperature_K=(t+1))
    dyn.run(steps_per_degree)
ase.io.write(os.path.join(end_folder, 'final.xyz'), atoms, append=False)
if equilibration_steps > 0:
    dyn.set_temperature(temperature_K=end_temp)
    dyn.run(equilibration_steps)
    ase.io.write(os.path.join(end_folder, 'final.xyz'), atoms, append=False)
print('MD finished!')

time_taken = time.time() - time_start
hours = int(time_taken / 3600)
time_taken -= hours * 3600
mins = int(time_taken / 60)
time_taken -= mins * 60
print(f'Wall time taken (hh:mm:ss) was {hours:02}:{mins:02}:{round(time_taken):02}.')

