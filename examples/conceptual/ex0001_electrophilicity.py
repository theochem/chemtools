
'''
================================
Chemical Potential & Hardness
================================

Compute chemical potential and hardness for CH4 using linear energy model.

'''

import os
from chemtools import ConceptualDFT_1File

# Temporary trick to find the data files
path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
file_path = os.path.join(path, 'data/test/ch4_uhf_ccpvdz.fchk')

# Build conceptual DFT descriptor tool
desp = ConceptualDFT_1File(file_path, model='linear')

# Print chemical potential and hardness
print 'Chemical Potential:', desp.mu
print 'Hardness          :', desp.eta
