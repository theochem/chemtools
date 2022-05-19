r"""
===================================================
EX1: Non-Covalent Interactions (NCI) of water dimer
===================================================

1. Build the NCI object for water dimmer using fchk files.
2. Dump files/scripts for visualizing NCI through VMD.
"""

from chemtools import NCI

# 1. Build NCI model using default settings

nci = NCI.from_file('../../data/test/h2o_dimer_pbe_sto3g.fchk')

# 2. Dump files/scripts for visualizing NCI
nci.dump_files('h2o_dimer')

# files generated are h2o_dimer-dens.cube, h2o_dimer-grad.cube, & h2o_dimer.vmd
# if VMD is setup on your system, you can visualize NCI with the command below, &
# obtain an image like the one displayed on the website above the script.
# $ vmd -e h2o_dimer.vmd