'''
===================================================
EX3: Non-Covalent Interactions (NCI) of water dimer
===================================================

The easiest way to calculate the Non-Covalent Interaction (NCI), by using the default settings is as follows:'''

from chemtools import NCI
nci = NCI.from_file('../../data/test/h2o_dimer_pbe_sto3g.fchk')
nci.dump_files('h2o_dimer')

# when you have vmd set up, you can visualize the files using "vmd -e h2o_dimer.vmd"
# this generates the output shown above.

# the code bolow is just to display the image
from tools.rug import plot_existing_image
plot_existing_image('nci_h2o_dimer.jpg')
