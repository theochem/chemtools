'''
==============================================
EX5: Plot Linear Fukui function (FMO Approach)
==============================================

1. Make a Cubic grid for plotting Fukui functions.
2. Build linear energy model for Formaldehyde, :math:`\mathbf{CH_2O}`,
   using frontier molecular orbital (FMO) theory.
3. Compute Fukui functions (f+, f- and f0) using linear energy model.
4. Make cube files & generate VMD (Visual Molecular Dynamics) scripts
   to plot Fukui function iso-surfaces.
'''

from chemtools import LocalConceptualDFT, CubeGen, print_vmd_script_isosurface, context

# 1. Make a Cubic grid for plotting Fukui functions.

# path to molecule's fchk file
file_path = context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk')
# make molecular cubic grid  with points spaced by 0.2 a.u. &
# extending 5.0 a.u. on every side of molecule
cube = CubeGen.from_file(file_path, spacing=0.2, threshold=5.0)

# 2. Build linear energy model for Formaldehyde using frontier molecular orbital (FMO) theory.

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
cube.dump_cube('coh2_ffm_fmo.cube', ffm)
cube.dump_cube('coh2_ffp_fmo.cube', ffp)
cube.dump_cube('coh2_ff0_fmo.cube', ff0)
# generate VMD scripts for visualizing iso-surfaces with VMD
print_vmd_script_isosurface('coh2_ffm_fmo.vmd', 'coh2_ffm_fmo.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2_ffp_fmo.vmd', 'coh2_ffp_fmo.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2_ff0_fmo.vmd', 'coh2_ff0_fmo.cube', isosurf=0.005)

# DISCARD BELOW:
# the code below is for displaying the ff image on the website, you should remove it
# when running the script on your machine.
from tools.rug import plot_existing_image
plot_existing_image('ch2o-ff_fmo.jpg')
