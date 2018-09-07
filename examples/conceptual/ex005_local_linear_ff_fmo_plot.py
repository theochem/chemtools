r"""
=========================================
EX5: Linear Fukui function (FMO Approach)
=========================================

1. Make a Cubic grid for plotting Fukui functions.
2. Build linear energy model for Formaldehyde, :math:`\mathbf{CH_2O}`,
   using frontier molecular orbital (FMO) theory.
3. Dump Fukui functions (f+, f- and f0) evaluated on cubic grid.
4. Generate VMD (Visual Molecular Dynamics) scripts to plot Fukui function iso-surfaces.
"""

from chemtools import LocalConceptualDFT, CubeGen, print_vmd_script_isosurface, context

file_path = context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk')

# 1. Make a Cubic grid for plotting Fukui functions.
#    The cubic grid points are spaced by 0.2 a.u. & extending 5.0 a.u. on each side.

cube = CubeGen.from_file(file_path, spacing=0.2, threshold=5.0)

# 2. Build linear energy model for Formaldehyde using frontier molecular orbital (FMO) theory.

tool = LocalConceptualDFT.from_file(file_path, model='linear', points=cube.points)

# 3. Dump Fukui functions (f+, f- and f0) evaluated on cubic grid.

cube.dump_cube('coh2_ffm_fmo.cube', tool.ff_minus)
cube.dump_cube('coh2_ffp_fmo.cube', tool.ff_plus)
cube.dump_cube('coh2_ff0_fmo.cube', tool.ff_zero)

# 4. Generate VMD scripts to plot Fukui function iso-surface.
#    To visualize the iso-surface, use command: $ vmd -e coh2_ffp_fmo.vmd

print_vmd_script_isosurface('coh2_ffp_fmo.vmd', 'coh2_ffp_fmo.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2_ff0_fmo.vmd', 'coh2_ff0_fmo.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2_ffm_fmo.vmd', 'coh2_ffm_fmo.cube', isosurf=0.005)

# DISCARD BELOW:
# the code below is for displaying the ff image on the website, you should remove it
# when running the script on your machine.
from tools.rug import plot_existing_image

plot_existing_image('ch2o-ff_fmo.jpg')
