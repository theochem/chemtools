'''
===================================================
EX5: Non-Covalent Interactions (NCI) of water dimer
===================================================

The easiest way to calculate the Non-Covalent Interaction (NCI), by using the default settings is as follows:'''

from chemtools import NCI
nci = NCI.from_file('../../data/test/h2o_dimer_pbe_sto3g.fchk')
nci.dump_files('h2o_dimer')

# when you have vmd set up, you can visualize the files using "vmd -e h2o_dimer.vmd"
# this generates the following output

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img=mpimg.imread('../../data/examples/images/nci_h2o_dimer.jpg')
imgplot = plt.imshow(img)
imgplot.axes.get_xaxis().set_visible(False)
imgplot.axes.get_yaxis().set_visible(False)

import os

os.system('rm *.vmd *.cube')
