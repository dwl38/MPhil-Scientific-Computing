#==================================================================================================
# MD simulation at 300 K, 1 atm over 500 ps
#==================================================================================================

import sys
import os
import time
import argparse
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
from torch.cuda import is_available as cuda_available

#==================================================================================================
# Modify parameters here

N_STEPS = 500000
TIMESTEP = 1 * units.fs

TTIME = 25 * units.fs
BULK_MODULUS = 2.2 * units.GPa
PTIME = 75 * units.fs
PFACTOR = PTIME * PTIME * BULK_MODULUS

#==================================================================================================
# Script

parser = argparse.ArgumentParser(prog='run.py')
parser.add_argument('--name', default=None)
parser.add_argument('--model', default=None)
parser.add_argument('--start', default='start.xyz')
parser.add_argument('--final', default='final.xyz')
args = parser.parse_args()

if args.name is None or args.model is None:
    print('[ERROR] Unspecified parameter')
    sys.exit()
os.mkdir(args.name)

DEVICE = 'cuda' if cuda_available() else 'cpu'
print(f'Loading MACE calculator "{args.model}" on "{DEVICE}"...')
calculator = MACECalculator(model_paths=args.model, device=DEVICE, enable_cueq=True)
atoms = ase.io.read(args.start)
atoms.calc = calculator

print('Starting MD...')
time_start = time.time()
MaxwellBoltzmannDistribution(atoms, temperature_K=300)
dyn = NPT(atoms, timestep=TIMESTEP, temperature_K=300, externalstress=1.01325 * units.bar,
          ttime=TTIME, pfactor=PFACTOR)
dyn.attach(MDLogger(dyn, atoms, os.path.join(args.name, 'md.log'), header=True, stress=True,
                    peratom=True, mode='w'), interval=100)
dyn.attach(ase.io.Trajectory(os.path.join(args.name, 'md.traj'), 'w', atoms).write, interval=100)

print('Running MD...')
for i in range(20):
    print(f'    Progress: {int(i * 5)}%...')
    dyn.run(int(N_STEPS/20))
    ase.io.write(os.path.join(args.name, args.final), atoms, append=False)
print('MD finished!')

time_taken = time.time() - time_start
hours = int(time_taken / 3600)
time_taken -= hours * 3600
mins = int(time_taken / 60)
time_taken -= mins * 60
print(f'Wall time taken (hh:mm:ss) was {hours:02}:{mins:02}:{round(time_taken):02}.')

