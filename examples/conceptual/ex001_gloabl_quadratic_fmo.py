'''
===========================================================
EX1: Global Quadratic Reactivity Descriptors (FMO Approach)
===========================================================

1. Build a quadratic energy model for :math:`\mathbf{CH_4}`
   using fontier molecular orbital (FMO) theory.
2. Print all global reactivity descriptors.
'''

from chemtools import ConceptualDFTGlobal

# 1. Build quadratic energy model

# relative path to molecule's file
file_path = '../../data/test/ch4_uhf_ccpvdz.fchk'
# build quadratic gloabl conceptual DFT tool (one file is passed, so FMO approach is taken)
tool = ConceptualDFTGlobal.from_file(file_path, model='quadratic')

# 2. Print all global reactivity descriptors

print 'Ionization Potential:', tool.ip
print 'Ionization Potential:', tool.ionization_potential
print 'Electron Affinity   :', tool.ea
print 'Electron Affinity   :', tool.electron_affinity
print 'Chemical Potential  :', tool.mu
print 'Chemical Potential  :', tool.chemical_potential
print 'Chemical Hardness   :', tool.eta
print 'Chemical Hardness   :', tool.chemical_hardness
print 'Chemical Softness   :', tool.softness
print 'Electronegativity   :', tool.electronegativity
print 'Electrophilicity    :', tool.electrophilicity
print 'Nucleophilicity     :', tool.nucleophilicity
print 'Nucleofugality      :', tool.nucleofugality
print 'N_max               :', tool.n_max
