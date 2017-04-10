r"""
=========================================================
EX10: 2D-Contours Quadratic Fukui Function (FMO Approach)
=========================================================

1. Make a 2D grid in the plane containing formaldehyde, :math:`\mathbf{CH_2O}`.
2. Build a quadratic energy models using frontier molecular orbital (FMO) theory approach.
3. Evaluate quadratic Fukui function on the grid points.
4. Plot 2D-contour plots of quadratic Fukui function.
"""

import numpy as np
import matplotlib.pyplot as plt
from chemtools import LocalConceptualDFT, context

# 1. Make a 2D grid in xy-plane (molecular plane).

# sample points in x and y direction (1D arrays); half-plane containing molecule is sampled
x = np.arange(-3., 4., 0.1)
y = np.arange(0., 3., 0.1)
z = np.zeros_like(x)

# make a mesh from 1D arrays
xv, yv, zv = np.meshgrid(x, y, z)

# turn mesh into 3D array
xyz = np.vstack((xv, yv, zv)).reshape((3, -1)).T
xyz = np.ascontiguousarray(xyz, dtype=np.float64)

# 2. Build a quadratic energy models using FMO approach

# path to molecule's fchk file
file_path = context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk')
# build a quadratic global conceptual DFT tool
tool = LocalConceptualDFT.from_file(file_path, model='quadratic', points=xyz)

# 3. Evaluate quadratic Fukui function on the grid points.

ff_quad = tool.fukui_function().reshape(xv.shape)[:, :, 0]

# 4. Plot 2D-contour plots of quadratic Fukui function.

# sample contour lines
levels = np.array([0.001 * n * n for n in range(500)])
# plot contours of quadratic Fukui function
plt.contour(x, y, ff_quad, levels)

# plot atomic centers
x_atoms, y_atoms = tool.coordinates[:, 0], tool.coordinates[:, 1]
plt.plot(x_atoms, y_atoms, marker='o', linestyle='', markersize=15, color='k')

# setting axis range & title
plt.axis([-3.0, 4.0, 0.0, 3.0])
plt.title('Quadratic Fukui Function Contour Plots', fontweight='bold')
# add axis label
plt.xlabel('x-axis', fontsize=12)
plt.ylabel('y-axis', fontsize=12)
# show plot
plt.show()
