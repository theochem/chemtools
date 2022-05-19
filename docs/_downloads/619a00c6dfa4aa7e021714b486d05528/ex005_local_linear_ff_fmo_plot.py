r"""
=================================
EX5: Fukui Function (Frontier MO)
=================================

Compute Fukui function on a cubic grid based on the linear energy model using
frontier molecular orbital (FMO) approach, and generate visualization scripts.

"""

from chemtools import LocalConceptualDFT, UniformGrid, print_vmd_script_isosurface

file_path = 'ch2o_q+0.fchk'

# 1. Make a Cubic grid for plotting Fukui functions.
#    The cubic grid points are spaced by 0.2 a.u. & extending 5.0 a.u. on each side.

cube = UniformGrid.from_file(file_path, spacing=0.2, extension=5.0)

# 2. Build linear energy model for Formaldehyde using frontier molecular orbital (FMO) theory.

tool = LocalConceptualDFT.from_file(file_path, model='linear', points=cube.points)

# 3. Dump Fukui functions (f+, f- and f0) evaluated on cubic grid.

cube.generate_cube('coh2_ffm_fmo.cube', tool.ff_minus)
cube.generate_cube('coh2_ffp_fmo.cube', tool.ff_plus)
cube.generate_cube('coh2_ff0_fmo.cube', tool.ff_zero)

# 4. Generate VMD scripts to plot Fukui function iso-surface.
#    To visualize the iso-surface, use command: $ vmd -e coh2_ffp_fmo.vmd

print_vmd_script_isosurface('coh2_ffp_fmo.vmd', 'coh2_ffp_fmo.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2_ff0_fmo.vmd', 'coh2_ff0_fmo.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2_ffm_fmo.vmd', 'coh2_ffm_fmo.cube', isosurf=0.005)