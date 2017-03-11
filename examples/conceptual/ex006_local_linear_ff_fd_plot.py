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

from chemtools import LocalConceptualDFT, CubeGen, print_vmd_script_isosurface

# 1. Make a Cubic grid for plotting Fukui functions.

# relative path to molecule's file
file_path = '../../data/examples/coh2_q+0_ub3lyp_6311g.fchk'
# make molecular cubic grid  with points spaced by 0.2 a.u. &
# extending 5.0 a.u. on every side of molecule
cube = CubeGen.from_file(file_path, spacing=0.2, threshold=5.0)


# 2. Build linear energy model for Formaldehyde using finite difference (FD) approach.

# make a list of 3 molecules' files used in finite difference approach.
# these molecules have the same geoemetry (they just differ in the number of electrons
# and multiplicity), so the cubic grid based on the first molecule made in section 1 works
# for all of them.
filenames = ['../../data/examples/coh2_q+0_ub3lyp_6311g.fchk',
             '../../data/examples/coh2_q+1_ub3lyp_6311g.fchk',
             '../../data/examples/coh2_q-1_ub3lyp_6311g.fchk']
tool = LocalConceptualDFT.from_file(filenames, model='linear', points=cube.points)

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
print_vmd_script_isosurface('coh2_ffm_fd.vmd', 'coh2_ffm_fd.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2_ffp_fd.vmd', 'coh2_ffp_fd.cube', isosurf=0.005)
print_vmd_script_isosurface('coh2_ff0_fd.vmd', 'coh2_ff0_fd.cube', isosurf=0.005)

# the code below is for displaying the ff image on the website, you should remove it
# when running the script on your machine.
# from tools.rug import plot_existing_image
# plot_existing_image('ch2o-ffp_fd.jpg')
