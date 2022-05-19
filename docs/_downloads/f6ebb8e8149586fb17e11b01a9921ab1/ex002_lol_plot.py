r"""
======================================================
EX2: LOL from wave-function file and user-defined cube
======================================================

Compute LOL and visualize it for formamide.

"""

from chemtools import LOL, UniformGrid

# 1. Build UniformGrid and LOL model

cub = UniformGrid.from_file('formamide_q+0.fchk', spacing=0.1, extension=2.0)
lol = LOL.from_file('formamide_q+0.fchk', grid=cub, trans='inverse_rational', trans_k=1, trans_a=1)

# 2. Generate cube file(s) and script for visualizing LOL
#    Files generated are chonh2-lol.cube & chonh2.vmd
#    To visualize the iso-surface, use command: $ vmd -e chonh2.vmd

lol.generate_scripts('chonh2', isosurf=0.55)