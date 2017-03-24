'''
======================================
EX8: Plot Energy Models (FMO Approach)
======================================

1. Build a linear and quadratic energy models for formaldehyde, :math:`\mathbf{CH_2O}`,
   using fontier molecular orbital (FMO) theory approach.
2. Compute energy values for various number of electrons.
3. Plot energy vs. number of electrons.
4. Plot data points used for modeling energy.
'''

import numpy as np
import matplotlib.pyplot as plt
from chemtools import GlobalConceptualDFT, context

# 1. Build linear and quadratic energy models using FMO approach

# path to molecule's fchk file
file_path = context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk')
# build linear & quadratic gloabl conceptual DFT tool (one file is passed, so FMO approach is taken)
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

# plot linear enery model
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
plt.title('Frontier Molecular Orbital Approach', fontsize=16, fontweight='bold')
# add legend & remove legend frame
plt.legend(frameon=False, fontsize=12)
# show plot
plt.tight_layout()
plt.show()
