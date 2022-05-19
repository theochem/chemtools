r"""
================================
EX1: ELF from wave-function file
================================

Compute ELF and visualize it for formaldehyde.

"""

from chemtools import ELF

# 1. Build ELF model

elf = ELF.from_file('ch2o_q+0.fchk', trans='rational', trans_k=2, trans_a=1)

# 2. Generate cube file(s) and script for visualizing ELF
#    Files generated are ch2o_q+0-elf.cube & ch2o_q+0.vmd
#    To visualize the iso-surface, use command: $ vmd -e ch2o_q+0.vmd

elf.generate_scripts('ch2o_q+0', isosurf=0.8)