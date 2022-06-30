"""
=========================================
EX2: ESP from Molecule & user-defied cube
=========================================

Compute ESP and visualize it on electron density iso-surface for 2,6-dichloropyridine.

"""

from chemtools import  UniformGrid, Molecule, print_vmd_script_isosurface

# 1. Build Molecule

fname = 'dichloropyridine26_q+0'

mol = Molecule.from_file(fname + '.fchk')
cub = UniformGrid.from_molecule(mol, spacing=1.0, extension=5.0)

# 2. Generate cube files: fname_esp.cube & fname_rho.cube

espname = fname + '_esp.cube'
rhoname = fname + '_rho.cube'

cub.generate_cube(rhoname, mol.compute_density(cub.points))
cub.generate_cube(espname, mol.compute_esp(cub.points))

# 3. Generate vmd script: fname.vmd
#    To visualize the iso-surface, use command: $ vmd -e fname.vmd

print_vmd_script_isosurface(fname + '.vmd', rhoname, colorfile=espname, isosurf=0.003,
                            scalemin=-0.02, scalemax=0.02)


# DISCARD BELOW:
# the code below is for displaying the ELF image on the website, you should remove it
# when running the script on your machine.
from tools.rug import plot_existing_image

plot_existing_image('esp_dichloropyridine26_q+0.png')
