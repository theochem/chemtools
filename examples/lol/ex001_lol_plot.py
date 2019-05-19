r"""
================================
EX1: LOL from wave-function file
================================

Compute LOL and visualize it for formaldehyde.

"""

from chemtools import LOL

# 1. Build LOL model

lol = LOL.from_file('ch2o_q+0.fchk', trans='inverse_rational', trans_k=1, trans_a=1)

# 2. Generate cube file(s) and script for visualizing LOL
#    Files generated are ch2o_q+0-lol.cube & ch2o_q+0.vmd
#    To visualize the iso-surface, use command: $ vmd -e ch2o_q+0.vmd

lol.generate_scripts('ch2o_q+0', isosurf=0.55)

# DISCARD BELOW:
# the code below is for displaying the LOL image on the website, you should remove it
# when running the script on your machine.
from tools.rug import plot_existing_image

plot_existing_image('lol055_ch2o_q+0.jpg')
