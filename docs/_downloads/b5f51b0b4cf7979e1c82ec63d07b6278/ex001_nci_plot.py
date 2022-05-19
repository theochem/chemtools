r"""
================================
EX1: NCI from wave-function file
================================

Compute NCI and visualize it for water dimer.

"""

from chemtools import NCI


# 1. Build NCI model

nci = NCI.from_file('h2o_dimer.fchk')

# 2. Generate plot, cube file(s) and script for visualizing NCI
#    Files generated are h2o_dimer-dens.cube, h2o_dimer-grad.cube, & h2o_dimer.vmd
#    To visualize the iso-surface, use command: $ vmd -e h2o_dimer.vmd

# nci.generate_plot('h2o_dimer')
nci.generate_scripts('h2o_dimer')