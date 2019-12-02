r"""
============================================
EX1: Average Local Ionization Potential (IP)
============================================

Compute average local IP and visualize it on electron density iso-surface for 2,6-dichloropyridine.

"""


from chemtools import  UniformGrid, Molecule, print_vmd_script_isosurface, DFTBasedTool

# 1. Build Orbital Theory model

fname = 'dichloropyridine26_q+0'

mol = Molecule.from_file(fname + '.fchk')
cub = UniformGrid.from_molecule(mol, spacing=1.0, extension=5.0)
orb = DFTBasedTool.from_molecule(mol, cub.points)

# 2. Generate cube files: fname_esp.cube & fname_rho.cube

lip_a, lip_b = orb.average_local_ionization_energy

lipname = fname + '_lip.cube'
rhoname = fname + '_rho.cube'

cub.generate_cube(rhoname, mol.compute_density(cub.points))
cub.generate_cube(lipname, lip_a + lip_b)

# 3. Generate vmd script: fname.vmd
#    To visualize the iso-surface, use command: $ vmd -e fname.vmd

print_vmd_script_isosurface(fname + '.vmd', rhoname, colorfile=lipname, isosurf=0.002,
                            scalemin=-1.0, scalemax=mol.mo.homo_energy[0])


# DISCARD BELOW:
# the code below is for displaying the ELF image on the website, you should remove it
# when running the script on your machine.
from tools.rug import plot_existing_image

plot_existing_image('local_ip.jpg')
