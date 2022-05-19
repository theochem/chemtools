r"""
=====================================
EX8: Plot Energy Models (FD Approach)
=====================================

Compute linear and quadratic energy models for various number of electrons
using finite Difference (FD) approach and plotting E vs. N.

"""

import numpy as np
import matplotlib.pyplot as plt
from chemtools import GlobalConceptualDFT

# 1. Build linear and quadratic energy models using FD approach

# make list of path to 3 molecule's fchk files used in finite difference approach
file_path = ['ch2o_q+0.fchk', 'ch2o_q+1.fchk', 'ch2o_q-1.fchk']
# build linear & quadratic global conceptual DFT tool (3 files are given, so FD approach is taken)
tool_lin = GlobalConceptualDFT.from_file(file_path, model='linear')
tool_qua = GlobalConceptualDFT.from_file(file_path, model='quadratic')

# 2. Compute energy values for various number of electrons.

# get reference number of electrons, n0, from either linear or quadratic models
n0 = tool_lin.n0
# sample number of electrons around n0
n_values = np.arange(n0 - 1.1, n0 + 1.1, 0.1)
# compute linear & quadratic energy values for sampled number of electrons
energy_lin = [tool_lin.energy(n) for n in n_values]
energy_qua = [tool_qua.energy(n) for n in n_values]

# 3. Plot energy vs. number of electrons.

# plot linear energy model
plt.plot(n_values, energy_lin, color='b', linestyle='-', linewidth=3,
         label='%s Model' % tool_lin.model.capitalize())
# plot quadratic energy model
plt.plot(n_values, energy_qua, color='r', linestyle='-', linewidth=3,
         label='%s Model' % tool_qua.model.capitalize())

# 4. Plot data points used for modeling energy.

# number of electrons used for modeling energy
n_data = [n0 - 1, n0, n0 + 1]
# compute energy values used for modeling energy
# Note: any of the tools built above can be used for this purpose because
#       they all have the same energy for N0 - 1, N0, and N0 + 1 electrons.
e_data = [tool_lin.energy(n) for n in n_data]

# plot given data points
plt.plot(n_data, e_data, marker='o', markersize=8, color='k', linestyle='', label='Given Data')

# add axis label
plt.xlabel('Number of electrons, N', fontsize=12)
plt.ylabel('Energy, E(N)', fontsize=12)
# add title
plt.title('Finite Difference Approach', fontsize=16, fontweight='bold')
# add legend & remove legend frame
plt.legend(frameon=False, fontsize=12)
# show plot
plt.tight_layout()
plt.show()