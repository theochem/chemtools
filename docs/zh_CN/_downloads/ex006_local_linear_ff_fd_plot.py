r"""
=============================================
EX6: Plot Linear Fukui function (FD Approach)
=============================================

1. Make a Cubic grid for plotting Fukui functions.
2. Build linear energy model for Formaldehyde, :math:`\mathbf{CH_2O}`,
   using finite difference (FD) approach.
3. Compute Fukui functions, i.e. f+, using linear energy model.
4. Make cube file & generate VMD (Visual Molecular Dynamics) script
   to plot Fukui function iso-surfaces.
"""

from chemtools import LocalConceptualDFT, CubeGen, print_vmd_script_isosurface, context

# 1. Make a Cubic grid for plotting Fukui functions.

# make list of path to 3 molecule's fchk files used in finite difference approach.
file_path = [context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk'),
             context.get_fn('examples/ch2o_q+1_ub3lyp_augccpvtz.fchk'),
             context.get_fn('examples/ch2o_q-1_ub3lyp_augccpvtz.fchk')]

# make molecular cubic grid  with points spaced by 0.2 a.u. &
# extending 5.0 a.u. on every side of molecule
# all 3 molecules have the same geometry (they just differ in the number of electrons
# and multiplicity), so the cubic grid based on the first molecule works for all.
cube = CubeGen.from_file(file_path[0], spacing=0.2, threshold=5.0)

# 2. Build linear energy model for Formaldehyde using finite difference (FD) approach.

# file_path contains 3 files are given, so FD approach is taken
tool = LocalConceptualDFT.from_file(file_path, model='linear', points=cube.points)

# 3. Compute Fukui functions (f+, f- and f0) using linear energy model.

# Fukui function minus
ffm = tool.ff_minus

# 4. Make cube files & generate VMD scripts to plot Fukui function iso-surfaces.

# dump Fukui function cubes file
cube.dump_cube('coh2_ffm_fd.cube', ffm)
# generate VMD scripts for visualizing iso-surfaces with VMD
print_vmd_script_isosurface('coh2_ffm_fd.vmd', 'coh2_ffm_fd.cube', isosurf=0.005, negative=True,
                            colorscheme=[0, 1])