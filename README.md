ChemTools
=========

[![Python](https://img.shields.io/badge/python-2.7-blue.svg)](https://docs.python.org/2.7/)
[![codecov](https://codecov.io/gh/QuantumElephant/chemtools/branch/master/graph/badge.svg?token=s2f4Ilawut)](https://codecov.io/gh/QuantumElephant/chemtools)
[![Build Status](https://travis-ci.com/QuantumElephant/chemtools.svg?token=wtCKs521Yw1urAV4F5DM&branch=master)](https://travis-ci.com/QuantumElephant/chemtools)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://github.com/QuantumElephant/chemtools/blob/master/LICENSE)


About
-----

<a href='https://chemtools.org'> ChemTools</a> is a free and open source Python library for interpreting the results of quantum
chemistry calculations. The goal of ChemTools is to provide a toolbox by which the
quantitative output of electronic structure theory calculations can be expressed in chemical
language. ChemTools provides easy-to-use core functionality to compute fundamental descriptors
of conceptual quantum chemistry, together with a flexible set of utilities allowing scientists
to easily test their own discoveries.

Website: https://chemtools.org


License
-------

ChemTools is distributed under GPL License version 3 (GPLv3)

Dependencies
------------

The following dependencies will be necessary for ChemTools to build properly,

* Python >= 2.7, < 3.0: http://www.python.org/ (Also install development files.)
* PIP >= 7.0: https://pip.pypa.io/ (Not required in some scenarios but never bad to have.)
* SciPy >= 0.11.0: http://www.scipy.org/
* NumPy >= 1.9.1: http://www.numpy.org/
* SymPy >= 0.7.1: http://code.google.com/p/sympy/
* Matplotlib >= 1.0: http://matplotlib.org/
* Nosetests >= 1.1.2: http://readthedocs.org/docs/nose/en/latest/
* HORTON >= 2.0.1: http://theochem.github.io/horton/2.0.1/index.html
* LFS >= 2.0.1: https://git-lfs.github.com/

Installation
------------

To install ChemTools :

```bash
python ./setup install --user
```

Check our website for more detailed
<a href='https://chemtools.org/usr_doc_installation.html'>installation guide</a>

Testing
-------

To run tests:

```bash
nosetests -v chemtools
```

Development
-----------
New contributors of all programming levels are welcome to join us. You can follow
our <a href='https://chemtools.org/tech_dev.html'>developer guidelines</a> for detailed information about contributing code, building
documentation and quality assurance.

Reference
---------
In anticipation of the first announcement of ChemTools in a scientific journal, please reference ChemTools as follows:
* F. Heidar-Zadeh, M. Richer, S. Fias, R.A. Miranda-Quintana, M. Chan,
M. Franco-Perez, C. E. Gonzalez-Espinoza, T.D. Kim, C. Lanssens,
A.H.G. Patel, X.D. Yang, E. Vohringer-Martinez, C. Cardenas, T. Verstraelen,
and P. W. Ayers. An explicit approach to conceptual density functional theory
descriptors of arbitrary order. Chem. Phys. Lett., 660:307–312, 2016. <a href='http://www.sciencedirect.com/science/article/pii/S0009261416305280'>
doi:10.1016/j.cplett.2016.07.039</a>.
