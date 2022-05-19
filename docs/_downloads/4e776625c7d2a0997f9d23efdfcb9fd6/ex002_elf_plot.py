r"""
======================================================
EX2: ELF from wave-function file and user-defined cube
======================================================

Compute ELF and visualize it for water dimer.

"""

from chemtools import ELF, UniformGrid

# 1. Build UniformGrid and ELF model

cub = UniformGrid.from_file('h2o_dimer.fchk', spacing=0.1, extension=2.0)
elf = ELF.from_file('h2o_dimer.fchk', grid=cub, trans='rational', trans_k=2, trans_a=1)

# 2. Generate cube file(s) and script for visualizing ELF
#    Files generated are h2o_dimer-elf.cube & h2o_dimer.vmd
#    To visualize the iso-surface, use command: $ vmd -e h2o_dimer.vmd

elf.generate_scripts('h2o_dimer', isosurf=0.8)