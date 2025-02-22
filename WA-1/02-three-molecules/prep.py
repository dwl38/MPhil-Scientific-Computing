import numpy as np
import ase.io
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

first_molecule = ase.io.read('../00-common/geoms/TATB.xyz')
second_molecule = first_molecule.copy()
second_molecule.translate(np.array((0.0, 0.0, 3.0)))
third_molecule = first_molecule.copy()
third_molecule.translate(np.array((0.0, 0.0, -3.0)))

atoms = first_molecule + second_molecule + third_molecule
atoms.set_cell(9.0 * np.identity(3))
atoms.set_pbc(True)

MaxwellBoltzmannDistribution(atoms, temperature_K=300)
ase.io.write('start.xyz', atoms)

