r"""
=============================================
EX7: Quadratic Dual Descriptor (FMO Approach)
=============================================

1. Make cubic grid for plotting dual descriptor.
2. Build quadratic energy model for Formaldehyde, :math:`\mathbf{CH_2O}`,
   using frontier molecular orbital (FMO) theory approach.
3. Dump dual descriptor evaluated on cubic grid.
4. Generate VMD (Visual Molecular Dynamics) script to visualize dual descriptor iso-surface.
"""

from chemtools import LocalConceptualDFT, CubeGen, print_vmd_script_isosurface, context

# 1. Make cubic grid for plotting dual descriptor.
#    The cubic grid points are spaced by 0.2 a.u. & extending 5.0 a.u. on each side.

file_path = context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk')
cube = CubeGen.from_file(file_path, spacing=0.2, threshold=5.0)

# 2. Build quadratic energy model for Formaldehyde using FMO approach.

tool = LocalConceptualDFT.from_file(file_path, model='quadratic', points=cube.points)

# 3. Dump dual descriptor evaluated on cubic grid.

cube.dump_cube('coh2_dual_fmo.cube', tool.dual_descriptor)

# 4. Generate VMD scripts to plot dual-descriptor iso-surface.
#    To visualize the iso-surface, use command: $ vmd -e coh2_dual_fmo.vmd

print_vmd_script_isosurface('coh2_dual_fmo.vmd', 'coh2_dual_fmo.cube', isosurf=0.005,
                            scalemin=-0.005, scalemax=0.005, colorscheme=[0, 1], negative=True)

# DISCARD BELOW:
# the code below is for displaying the dual descriptor image on the website, you should remove it
# when running the script on your machine.
from tools.rug import plot_existing_image

plot_existing_image('ch2o_dual_fmo.jpg')
