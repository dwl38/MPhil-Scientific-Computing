#==================================================================================================
# MD simulation of gradual heating and then thermal equilibration at a user-specified temperature,
# under constant pressure of 1 atm, for measurements of the lattice parameters and thermal
# expansion coefficients.
#==================================================================================================

import sys
import os
import shutil
import time
import numpy as np
import ase.io

# To suppress warnings from PyTorch coming in through MACE
if not sys.warnoptions:
    import warnings
    warnings.simplefilter('ignore')

from mace.calculators import MACECalculator
from ase.md import MDLogger
from ase.md.npt import NPT
from ase import units
from torch.cuda import is_available as cuda_available

#==================================================================================================
# Modify parameters here

TIMESTEP = 1 * units.fs

PATM = 1.01325 * units.bar

TTIME = 25 * units.fs
BULK_MODULUS = 100 * units.GPa
PTIME = 75 * units.fs
PFACTOR = PTIME * PTIME * BULK_MODULUS

MACE_MODEL = '../00-common/mace-models/MACE-OFF24_medium.model'

#==================================================================================================
# Script

start_temp = float(sys.argv[1]) if len(sys.argv) > 1 else 300
end_temp = float(sys.argv[2]) if len(sys.argv) > 2 else 400
temp_incre = float(sys.argv[3]) if len(sys.argv) > 3 else 1
steps_per_degree = int(sys.argv[4]) if len(sys.argv) > 4 else 0
equilibration_steps = int(sys.argv[5]) if len(sys.argv) > 5 else 0
time_start = time.time()
print(f'[{time.ctime()}]: Running crystal gradual heat from {start_temp}K to {end_temp}K with' +
      f' {steps_per_degree} timesteps per {temp_incre} degrees, and then {equilibration_steps}' +
      ' final equilibration timesteps...')

start_folder = str(int(round(start_temp))) + 'K'
end_folder = str(int(round(end_temp))) + 'K'
if not os.path.isdir(start_folder) or not os.path.isfile(os.path.join(start_folder, 'final.xyz')):
    print(f'[ERROR] Start folder {start_folder} does not exist!')
    sys.exit()
if os.path.isdir(end_folder):
    shutil.rmtree(end_folder)
    print(f'[WARNING] Pre-existing end folder "{end_folder}" deleted!')
os.mkdir(end_folder)

DEVICE = 'cuda' if cuda_available() else 'cpu'
print(f'Loading MACE calculator on "{DEVICE}"...')
calculator = MACECalculator(model_paths=MACE_MODEL, device=DEVICE, enable_cueq=True)
atoms = ase.io.read(os.path.join(start_folder, 'final.xyz'))
atoms.calc = calculator

if steps_per_degree > 0:

    print('Running MD for gradual heating process...')
    dyn = NPT(atoms, timestep=TIMESTEP, temperature_K=start_temp, externalstress=PATM,
              ttime=TTIME, pfactor=PFACTOR)
    dyn.attach(MDLogger(dyn, atoms, os.path.join(end_folder, 'heat.log'),
               header=True, stress=True, peratom=True, mode='w'), interval=100)
    dyn.attach(ase.io.Trajectory(os.path.join(end_folder, 'heat.traj'), 'w', atoms).write,
               interval=10)
    n_steps = int((end_temp - start_temp) / temp_incre) + 1
    for t in np.linspace(start_temp, end_temp, n_steps):
        print(f'    Progress: running {t}K...')
        dyn.set_temperature(temperature_K=t)
        dyn.run(steps_per_degree)
    ase.io.write(os.path.join(end_folder, 'final.xyz'), atoms, append=False)

if equilibration_steps > 0:

    print('Running MD for equilibration...')
    dyn = NPT(atoms, timestep=TIMESTEP, temperature_K=end_temp, externalstress=PATM,
              ttime=TTIME, pfactor=PFACTOR)
    dyn.attach(MDLogger(dyn, atoms, os.path.join(end_folder, 'equi.log'),
               header=True, stress=True, peratom=True, mode='w'), interval=100)
    dyn.attach(ase.io.Trajectory(os.path.join(end_folder, 'equi.traj'), 'w', atoms).write,
               interval=10)
    for i in range(10):
        print(f'    Progress: {i * 10}%...')
        dyn.run(int(equilibration_steps / 10))
    ase.io.write(os.path.join(end_folder, 'final.xyz'), atoms, append=False)
print('MD finished!')

time_taken = time.time() - time_start
hours = int(time_taken / 3600)
time_taken -= hours * 3600
mins = int(time_taken / 60)
time_taken -= mins * 60
print(f'[{time.ctime()}]: Wall time taken was {hours:02}:{mins:02}:{round(time_taken):02}.')

