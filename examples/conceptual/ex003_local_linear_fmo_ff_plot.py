'''
==============================================
EX3: Plot Linear Fukui function (FMO Approach)
==============================================

1. Make a Cubic grid for plotting Fukui functions.
2. Build linear energy model for Formaldehyde, :math:`\mathbf{CH_2O}`,
   using frontier molecular orbital (FMO) theory.
3. Compute Fukui functions (f+, f- and f0) using linear energy model.
4. Generate VMD (Visual Molecular Dynamics) scripts to plot Fukui function
   iso-surfaces.
'''

from chemtools import LocalConceptualDFT, CubeGen, print_vmd_script_isosurface

# 1. Make a Cubic grid for plotting Fukui functions.

# relative path to molecule's file
file_path = '../../data/examples/coh2_q+0_ub3lyp_6311g.fchk'
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

# 4. Generate VMD (Visual Molecular Dynamics) scripts to plot Fukui function iso-surfaces.

# dump Fukui function cubes files
cube.dump_cube('coh2-ffm.cube', ffm)
cube.dump_cube('coh2-ffp.cube', ffp)
cube.dump_cube('coh2-ff0.cube', ff0)
# construct VMD scripts for visualizing the iso-surfaces with VMD
print_vmd_script_isosurface('coh2-ffm.vmd', 'coh2-ffm.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2-ffp.vmd', 'coh2-ffp.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2-ff0.vmd', 'coh2-ff0.cube', isosurf=0.005)
