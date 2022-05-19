r"""
============================================
EX3: ELF from Molecule and user-defined cube
============================================

Compute ELF and visualize it for formamide.

"""

from chemtools import ELF, UniformGrid, Molecule

# 1. Build Molecule, UnifromGrid and ELF model

mol = Molecule.from_file('formamide_q+0.fchk')
cub = UniformGrid.from_molecule(mol, spacing=0.1, extension=2.0)
elf = ELF.from_molecule(mol, grid=cub, trans='hyperbolic', trans_k=1, trans_a=1)

# 2. Generate cube file(s) and script for visualizing ELF
#    Files generated are chonh2-elf.cube & chonh2.vmd
#    To visualize the iso-surface, use command: $ vmd -e chonh2.vmd

elf.generate_scripts('chonh2', isosurf=0.8)