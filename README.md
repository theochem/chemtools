ChemTools
=========
<a href='https://travis-ci.com/QuantumElephant/chemtools'><img src='https://travis-ci.com/QuantumElephant/chemtools.svg?token=wtCKs521Yw1urAV4F5DM&branch=master'></a>
<a href='https://docs.python.org/2.7/'><img src='https://img.shields.io/badge/python-2.7-blue.svg'></a>

Refer <a href='http://chemtools.org'> ChemTools</a> for Detailed Documentation and installation Guide

About
-----

<a href='http://chemtools.org'> ChemTools</a> is a free and open source Python library for interpreting the results of quantum
chemistry calculations. The goal of ChemTools is to provide a toolbox by which the
quantitative output of electronic structure theory calculations can be expressed in chemical
language. ChemTools provides easy-to-use core functionality to compute fundamental descriptors
of conceptual quantum chemistry, together with a flexible set of utilities allowing scientists
to easily test their own discoveries.

License
-------

ChemTools is distributed under GPL License version 3 (GPLv3)

Dependence
----------

The following dependencies will be necessary for ChemTools to build properly,

* Python >= 2.7, < 3.0: http://www.python.org/ (Also install development files.)
* PIP >= 7.0: https://pip.pypa.io/ (Not required in some scenarios but never bad to have.)
* SciPy >= 0.11.0: http://www.scipy.org/
* NumPy >= 1.9.1: http://www.numpy.org/
* Sympy >= 0.7.1: http://code.google.com/p/sympy/
* Matplotlib >= 1.0: http://matplotlib.org/
* Nosetests >= 1.1.2: http://readthedocs.org/docs/nose/en/latest/
* HORTON >= 2.0.1: http://theochem.github.io/horton/2.0.1/index.html
* LFS: https://git-lfs.github.com/

Installation
------------

To install ChemTools :

```bash
python ./setup install --user
```

Tests
-----

To run tests:

```bash
nosetests -v saddle
```
