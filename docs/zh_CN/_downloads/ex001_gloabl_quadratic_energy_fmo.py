r"""
===========================================================
EX1: Global Quadratic Reactivity Descriptors (FMO Approach)
===========================================================

1. Build a quadratic energy model for formaldehyde, :math:`\mathbf{CH_2O}`,
   using frontier molecular orbital (FMO) theory approach.
2. Print all available global quadratic reactivity descriptors.
3. Compute quadratic energy model and its derivatives for various numbers of electrons.
"""

from chemtools import GlobalConceptualDFT, context

# 1. Build quadratic energy model using FMO approach

# path to molecule's fchk file
file_path = context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk')
# build quadratic global conceptual DFT tool (one file is passed, so FMO approach is taken)
tool = GlobalConceptualDFT.from_file(file_path, model='quadratic')

# 2. Print all available global quadratic reactivity descriptors

print 'Ionization Potential:', tool.ip, tool.ionization_potential
print 'Electron Affinity   :', tool.ea, tool.electron_affinity
print 'Chemical Potential  :', tool.mu, tool.chemical_potential
print 'Chemical Hardness   :', tool.eta, tool.chemical_hardness
print 'Chemical Softness   :', tool.softness
print 'Electronegativity   :', tool.electronegativity
print 'Electrophilicity    :', tool.electrophilicity
print 'Electrofugality     :', tool.electrofugality
print 'Nucleofugality      :', tool.nucleofugality
print 'N_max               :', tool.n_max
print

# 3. Compute quadratic energy model and its derivatives for various numbers of electrons.

print 'Energy at N=15.5:', tool.energy(15.5)
print 'Energy at N=16.0:', tool.energy(16.0)
print 'Energy at N=16.5:', tool.energy(16.5)
print 'Energy at N=Nmax:', tool.energy(tool.n_max)
print
print '1st Derivative Energy at N=15.5:', tool.energy_derivative(15.5, order=1)
print '1st Derivative Energy at N=16.0:', tool.energy_derivative(16.0, order=1)
print '1st Derivative Energy at N=16.5:', tool.energy_derivative(16.5, order=1)
print '1st Derivative Energy at N=Nmax:', tool.energy_derivative(tool.n_max, order=1)
print
print '2nd Derivative Energy at N=15.5:', tool.energy_derivative(15.5, order=2)
print '2nd Derivative Energy at N=16.0:', tool.energy_derivative(16.0, order=2)
print '2nd Derivative Energy at N=16.5:', tool.energy_derivative(16.5, order=2)
print '2nd Derivative Energy at N=Nmax:', tool.energy_derivative(tool.n_max, order=2)
print