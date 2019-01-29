r"""
===================================================
EX1: Non-Covalent Interactions (NCI) of water dimer
===================================================
"""

from chemtools import NCI


# 1. Build NCI model using default settings

nci = NCI.from_file('data/examples/h2o_dimer_pbe_sto3g.fchk')

# 2. Dump files/scripts for visualizing NCI
#    Files generated are h2o_dimer-dens.cube, h2o_dimer-grad.cube, & h2o_dimer.vmd
#    To visualize the iso-surface, use command: $ vmd -e h2o_dimer.vmd

nci.dump_files('h2o_dimer')


# DISCARD BELOW:
# the code below is for displaying the NCI image on the website, you should remove it
# when running the script on your machine.
from tools.rug import plot_existing_image

plot_existing_image('nci_h2o_dimer.jpg')
