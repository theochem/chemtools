ChemTools
=========

[![codecov](https://codecov.io/gh/QuantumElephant/chemtools/branch/master/graph/badge.svg?token=s2f4Ilawut)](https://codecov.io/gh/QuantumElephant/chemtools)
[![Build Status](https://travis-ci.com/QuantumElephant/chemtools.svg?token=wtCKs521Yw1urAV4F5DM&branch=master)](https://travis-ci.com/QuantumElephant/chemtools)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://github.com/QuantumElephant/chemtools/blob/master/LICENSE)


About
-----
[ChemTools](https://chemtools.org)  is a free and open source Python library for interpreting the results of quantum
chemistry calculations. The goal of ChemTools is to provide a toolbox by which the
quantitative output of electronic structure theory calculations can be expressed in chemical
language. ChemTools provides easy-to-use core functionality to compute fundamental descriptors
of conceptual quantum chemistry, together with a flexible set of utilities allowing scientists
to easily test their own discoveries. Visit [ChemTools website](https://chemtools.org) for more information.

Citation
--------
Please use the following citation when using ChemTools in your research:

> L. Pujal, A. Tehrani, and F. Heidar-Zadeh. ChemTools: Gain Chemical
  Insight from Quantum Chemistry Calculations. In Conceptual Density Functional Theory: Towards A New Chemical Reactivity Theory (Editor: Shubin Liu), [Wiley, 2022; pp 649-661](https://onlinelibrary.wiley.com/doi/abs/10.1002/9783527829941.ch32).

> F. Heidar-Zadeh, M. Richer, S. Fias, R.A. Miranda-Quintana, M. Chan,
  M. Franco-Perez, C. E. Gonzalez-Espinoza, T.D. Kim, C. Lanssens, A.H.G. Patel, X.D. Yang, E. Vohringer-Martinez, C. Cardenas, T. Verstraelen, and P. W. Ayers. An explicit approach to conceptual density functional theory descriptors of arbitrary order. [Chem. Phys. Lett., 660:307–312, 2016.](http://www.sciencedirect.com/science/article/pii/S0009261416305280)

Installation
------------
To install ChemTools and its dependencies using Anaconda,
```bash
# make a Python environment (with any tool you use, pyenv in this guide) and activate it
pyenv virtualenv chemtools_py313
pyenv activate chemtools_py313

# install HORTON 3 libraries
pip install git+https://github.com/theochem/iodata.git
pip install git+https://github.com/theochem/gbasis.git
pip install git+https://github.com/theochem/grid.git
pip install git+https://github.com/theochem/denspart.git

# install the correct version of dependencies
pip install numpy scipy==1.16.2 matplotlib

# clone the repository
git clone https://github.com/theochem/chemtools.git

# install Chemtools
cd chemtools
pip install -e .
```

Check our website for more detailed [installation guide](https://chemtools.org/usr_doc_installation.html).

Development
-----------
New contributors of all programming levels are welcome to join us. You can follow our [developer guidelines](https://chemtools.org/tech_dev.html) for detailed information about contributing code, building
documentation and quality assurance.
