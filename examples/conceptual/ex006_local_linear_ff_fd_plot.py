'''
=============================================
EX6: Plot Linear Fukui function (FD Approach)
=============================================

1. Make a Cubic grid for plotting Fukui functions.
2. Build linear energy model for Formaldehyde, :math:`\mathbf{CH_2O}`,
   using finite difference (FD) approach.
3. Compute Fukui functions (f+, f- and f0) using linear energy model.
4. Make cube files & generate VMD (Visual Molecular Dynamics) scripts
   to plot Fukui function iso-surfaces.
'''

from chemtools import LocalConceptualDFT, CubeGen, print_vmd_script_isosurface, context

# 1. Make a Cubic grid for plotting Fukui functions.

# make list of path to 3 molecule's fchk files used in finite difference approach.
file_path = [context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk'),
             context.get_fn('examples/ch2o_q+1_ub3lyp_augccpvtz.fchk'),
             context.get_fn('examples/ch2o_q-1_ub3lyp_augccpvtz.fchk')]

# make molecular cubic grid  with points spaced by 0.2 a.u. &
# extending 5.0 a.u. on every side of molecule
# all 3 molecules have the same geoemetry (they just differ in the number of electrons
# and multiplicity), so the cubic grid based on the first molecule works for all.
cube = CubeGen.from_file(file_path[0], spacing=0.2, threshold=5.0)

# 2. Build linear energy model for Formaldehyde using finite difference (FD) approach.

# file_path contains 3 files are given, so FD approach is takn
tool = LocalConceptualDFT.from_file(file_path, model='linear', points=cube.points)

# 3. Compute Fukui functions (f+, f- and f0) using linear energy model.

# Fukui function minus
ffm = tool.ff_minus
# Fukui function plus
ffp = tool.ff_plus
# Fukui function zero
ff0 = tool.ff_zero

# 4. Make cube files & generate VMD scripts to plot Fukui function iso-surfaces.

# dump Fukui function cubes files
cube.dump_cube('coh2_ffm_fd.cube', ffm)
cube.dump_cube('coh2_ffp_fd.cube', ffp)
cube.dump_cube('coh2_ff0_fd.cube', ff0)
# generate VMD scripts for visualizing iso-surfaces with VMD
print_vmd_script_isosurface('coh2_ffm_fd.vmd', 'coh2_ffm_fd.cube', isosurf=0.005, negative=True, colorscheme=[0, 1])
print_vmd_script_isosurface('coh2_ffp_fd.vmd', 'coh2_ffp_fd.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2_ff0_fd.vmd', 'coh2_ff0_fd.cube', isosurf=0.005)

# DISCARD BELOW:
# the code below is for displaying the ff image on the website, you should remove it
# when running the script on your machine.
from tools.rug import plot_existing_image
plot_existing_image('ch2o-ffm_fd.jpg')
