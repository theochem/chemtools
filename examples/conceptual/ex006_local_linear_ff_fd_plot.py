r"""
========================================
EX6: Linear Fukui function (FD Approach)
========================================

1. Make a Cubic grid for plotting Fukui functions.
2. Build linear energy model for Formaldehyde, :math:`\mathbf{CH_2O}`,
   using finite difference (FD) approach.
3. Dump Fukui functions (f+, f- and f0) evaluated on cubic grid.
4. Generate VMD (Visual Molecular Dynamics) scripts to plot Fukui function iso-surfaces.
"""

from chemtools import LocalConceptualDFT, CubeGen, print_vmd_script_isosurface, context

file_path = [context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk'),
             context.get_fn('examples/ch2o_q+1_ub3lyp_augccpvtz.fchk'),
             context.get_fn('examples/ch2o_q-1_ub3lyp_augccpvtz.fchk')]

# 1. Make cubic grid for plotting Fukui function.
#    The cubic grid points are spaced by 0.2 a.u. & extending 5.0 a.u. on each side.
#    All 3 molecules have the same geometry (they just differ in the number of electrons
#    and multiplicity), so the cubic grid based on the first molecule works for all.

cube = CubeGen.from_file(file_path[0], spacing=0.2, threshold=5.0)

# 2. Build linear energy model for Formaldehyde using finite difference (FD) approach.

tool = LocalConceptualDFT.from_file(file_path, model='linear', points=cube.points)

# 3. Dump Fukui functions (f+, f- and f0) evaluated on cubic grid.

cube.dump_cube('coh2_ffp_fd.cube', tool.ff_plus)
cube.dump_cube('coh2_ff0_fd.cube', tool.ff_zero)
cube.dump_cube('coh2_ffm_fd.cube', tool.ff_minus)

# 4. Generate VMD scripts to plot dual-descriptor iso-surface.
#    To visualize the iso-surface, use command: $ vmd -e coh2_ffp_fd.vmd

print_vmd_script_isosurface('coh2_ffp_fd.vmd', 'coh2_ffp_fd.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2_ff0_fd.vmd', 'coh2_ff0_fd.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2_ffm_fd.vmd', 'coh2_ffm_fd.cube', isosurf=0.005,
                            negative=True, colorscheme=[0, 1])

# DISCARD BELOW:
# the code below is for displaying the ff image on the website, you should remove it
# when running the script on your machine.
from tools.rug import plot_existing_image

plot_existing_image('ch2o-ffm_fd.jpg')
