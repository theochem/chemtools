'''
==================================================
EX7: Plot Quadratic Dual Descriptor (FMO Approach)
==================================================

1. Make cubic grid for plotting dual descriptor.
2. Build quadratic energy model for Formaldehyde, :math:`\mathbf{CH_2O}`,
   using frontier molecular orbitral (FMO) theory approach.
3. Compute dual descriptor using quadratic energy model.
4. Make dual descriptor cube file & generate VMD (Visual Molecular Dynamics)
   script to visualize its iso-surface.
'''

from chemtools import LocalConceptualDFT, CubeGen, print_vmd_script_isosurface

# 1. Make cubic grid for plotting dual descriptor.

# relative path to molecule's file
file_path = '../../data/examples/coh2_q+0_ub3lyp_6311g.fchk'
# make molecular cubic grid  with points spaced by 0.2 a.u. &
# extending 5.0 a.u. on every side of molecule
cube = CubeGen.from_file(file_path, spacing=0.2, threshold=5.0)

# 2. Build quadratic energy model for Formaldehyde using FMO approach.

tool = LocalConceptualDFT.from_file(file_path, model='quadratic', points=cube.points)

# 3. Compute dual descriptor using quadratic energy model.

dual = tool.dual_descriptor()

# 4. Make dual descriptor cube file & generate VMD scripts to plot its iso-surface.

# dump dual descriptor cube file
cube.dump_cube('coh2_dual_fmo.cube', dual)
# generate VMD scripts for visualizing iso-surfaces with VMD
print_vmd_script_isosurface('coh2_dual_fmo.vmd', 'coh2_dual_fmo.cube', isosurf=0.005,
                            scalemin=-0.005, scalemax=0.005, colorscheme=[0, 1], negative=True)

# the code below is for displaying the dual descriptor image on the website, you should remove it
# when running the script on your machine.
from tools.rug import plot_existing_image
plot_existing_image('ch2o_dual_fmo.jpg')
