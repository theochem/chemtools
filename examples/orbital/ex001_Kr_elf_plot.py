'''
======================================================
EX1: Plot Krypton Electron Localization Function (ELF)
======================================================

1. Run a Hartree-Fock calculation with cc-pVDZ basis-set using HORTON.
2. Build an OrbitalLocalTool.
3. Plot Electron Localization Function (ELF) & radially weighted density.
'''

import numpy as np
import matplotlib.pyplot as plt
from horton import *
from chemtools import OrbitalLocalTool

# 1. Run a Hartree-Fock calculation with cc-pVDZ basis-set using HORTON

# construct atom object placed at origin
mol = IOData(title='Krypton')
mol.coordinates = np.array([[0.0, 0.0, 0.0]])
mol.numbers = np.array([36])

# make a cc-pVDZ basis-set
obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pVDZ')

# creata linalg factory
lf = CholeskyLinalgFactory(obasis.nbasis)

# compute Gaussian integrals
olp = obasis.compute_overlap(lf)
kin = obasis.compute_kinetic(lf)
na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
er = obasis.compute_electron_repulsion(lf)

# create alpha orbitals
exp_alpha = lf.create_expansion()

# initial guess
guess_core_hamiltonian(olp, kin, na, exp_alpha)

# construct restricted HF effective Hamiltonian
external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
terms = [RTwoIndexTerm(kin, 'kin'),
         RDirectTerm(er, 'hartree'),
         RExchangeTerm(er, 'x_hf'),
         RTwoIndexTerm(na, 'ne'),
        ]
ham = REffHam(terms, external)

# select orbital occupation scheme (1 alpha electron)
occ_model = AufbauOccModel(18)

# compute WFN with plain SCF
scf_solver = PlainSCFSolver(1e-6)
scf_solver(ham, lf, olp, occ_model, exp_alpha)

# 2. Post-processing
# ------------------

# construct an atomic grid
grid = AtomicGrid(mol.numbers[0], mol.pseudo_numbers[0], mol.coordinates[0],
                  agspec='exp:0.005:10.0:200:194', random_rotate=False)

# get radial distribution
radii = grid.rgrid.radii

# construce an OrbitalLocalTool
tool = OrbitalLocalTool(grid.points, obasis, exp_alpha)

# calculate spherically-averaged density & Electron Localization Function (ELF)
dens = grid.get_spherical_average(tool.density)
elf = grid.get_spherical_average(tool.elf)

# calculate radially weighted density, i.e. $4.0 \pi r^2 \rho$
rdens = 4.0 * np.pi * radii * radii * dens

# 3. Plot Electron Localization Function (ELF) & radially weighted density.

# compute log10 of radii
logx = np.array([np.math.log(radius, 10.) for radius in grid.rgrid.radii])

# plot ElF & radially weighted density
plt.plot(logx, elf, linewidth=2, label='Electron Localization Function')
plt.plot(logx, 0.01875 * rdens, linewidth=2, label='Scaled radially weighted density')
plt.legend(loc='upper right', fontsize=10)
plt.xlim(-2., 1.)
plt.xlabel('log$_{10}$(radial distance from the nucleus / $a_0$ )')
plt.show()
