r"""
=================================
EX6: Fukui function (Finite Diff)
=================================

Compute Fukui function on a cubic grid based on the linear energy model using
finite difference (FD) approach, and generate visualization scripts.

"""

from chemtools import LocalConceptualDFT, UniformGrid, print_vmd_script_isosurface

file_path = ['ch2o_q+0.fchk', 'ch2o_q+1.fchk', 'ch2o_q-1.fchk']

# 1. Make cubic grid for plotting Fukui function.
#    The cubic grid points are spaced by 0.2 a.u. & extending 5.0 a.u. on each side.
#    All 3 molecules have the same geometry (they just differ in the number of electrons
#    and multiplicity), so the cubic grid based on the first molecule works for all.

cube = UniformGrid.from_file(file_path[0], spacing=0.2, extension=5.0)

# 2. Build linear energy model for Formaldehyde using finite difference (FD) approach.

tool = LocalConceptualDFT.from_file(file_path, model='linear', points=cube.points)

# 3. Dump Fukui functions (f+, f- and f0) evaluated on cubic grid.

cube.generate_cube('coh2_ffp_fd.cube', tool.ff_plus)
cube.generate_cube('coh2_ff0_fd.cube', tool.ff_zero)
cube.generate_cube('coh2_ffm_fd.cube', tool.ff_minus)

# 4. Generate VMD scripts to plot dual-descriptor iso-surface.
#    To visualize the iso-surface, use command: $ vmd -e coh2_ffp_fd.vmd

print_vmd_script_isosurface('coh2_ffp_fd.vmd', 'coh2_ffp_fd.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2_ff0_fd.vmd', 'coh2_ff0_fd.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2_ffm_fd.vmd', 'coh2_ffm_fd.cube', isosurf=0.005,
                            negative=True, colorscheme=[0, 1])