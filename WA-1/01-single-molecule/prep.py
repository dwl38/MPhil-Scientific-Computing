import numpy as np
import ase.io
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

atoms = ase.io.read('../00-common/geoms/TATB.xyz')

width_x = np.max(atoms.positions[:,0]) - np.min(atoms.positions[:,0])
width_y = np.max(atoms.positions[:,1]) - np.min(atoms.positions[:,1])
width_z = np.max(atoms.positions[:,2]) - np.min(atoms.positions[:,2])
width = 2 * max(max(width_x, width_y), width_z)

atoms.set_cell(width * np.identity(3))
atoms.set_pbc(False)

MaxwellBoltzmannDistribution(atoms, temperature_K=300)
ase.io.write('start.xyz', atoms)

