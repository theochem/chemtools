'''
===============================================
EX2: Plot Quadratic Energy Model (FMO Approach)
===============================================

1. Build a quadratic energy model for :math:`\mathbf{CH_4}` using fontier molecular orbital (FMO) theory.
2. Compute quadratic energy model for various number of electrons.
3. Plot energy versus number of electrons.
'''

import numpy as np
import matplotlib.pyplot as plt
from chemtools import GlobalConceptualDFT

# 1. Build quadratic energy model

# relative path to molecule's file
file_path = '../../data/test/ch4_uhf_ccpvdz.fchk'
# build quadratic gloabl conceptual DFT tool
tool = GlobalConceptualDFT.from_file(file_path, model='quadratic')

# 2. Compute energy for sampled number of electrons

# sample values for number of electrons, N, between 5.0 and 15.0 with step=0.5
n_values = np.arange(tool.n0 - 5.0, tool.n0 + 5.0, 0.5)
# compute energy values for the sampled number of electrons
e_values = [tool.energy(n) for n in n_values]

# 3. Plot energy versus number of electrons

# plot sampled points as a dashed line
plt.plot(n_values, e_values, linestyle='--', linewidth=3, color='b', label=tool.model)
# plot original points used for fitting energy model
n_values = [tool.n0 - 1, tool.n0, tool.n0 + 1]
e_values = [tool.energy(n) for n in n_values]
plt.plot(n_values, e_values, marker='o', markersize=6, color='r', linestyle='', label='Data')
# add axis label
plt.xlabel('Number of electrons, N')
plt.ylabel('Energy, E(N)')
# add legend & remove legend frame
plt.legend(frameon=False, fontsize=12)
plt.show()
