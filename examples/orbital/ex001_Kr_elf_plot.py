r"""
=================================================
EX1: Krypton Electron Localization Function (ELF)
=================================================
"""


import numpy as np
import matplotlib.pyplot as plt

from chemtools import Molecule, ELF
from horton import AtomicGrid


# load atom and construct and atomic grid
filename = 'data/examples/atom_kr_hf_ccpvdz.fchk'
molecule = Molecule.from_file(filename)
grid = AtomicGrid(molecule.numbers[0], molecule.pseudo_numbers[0], molecule.coordinates[0],
                  agspec='exp:0.005:10.0:200:194', random_rotate=False)

# construct an instance of ELF
tool = ELF.from_molecule(molecule, spin='ab', index=None, grid=grid)

# calculate spherically-averaged density ($4.0 \pi r^2 \rho$) & ELF
elf = grid.get_spherical_average(tool.values)
dens = grid.get_spherical_average(tool.density)
dens *= (4.0 * np.pi * grid.rgrid.radii**2)

# compute log10 of radii & plot ElF and (scaled) radially weighted density
logx = np.log10(grid.rgrid.radii)
plt.plot(logx, elf, linewidth=2, label='Electron Localization Function')
plt.plot(logx, 0.01875 * dens, linestyle='--', linewidth=2, label='Radially Weighted Density (Scaled)')
plt.legend(loc='upper right', fontsize=8, frameon=False)
plt.xlim(-2., 1.)
plt.xlabel('log$_{10}$(radial distance from the nucleus / $a_0$)')
plt.show()
