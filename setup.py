#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2014-2015 The ChemTools Development Team
#
# This file is part of ChemTools.
#
# ChemTools is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# ChemTools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --

from distutils.core import setup


setup(
    name='cdft',
    version='0.0',
    description='Conceptual Density Functional Theory (DFT) Reactivity Descriptors Package.',
    author='Farnaz Heidar-Zadeh',
    author_email='heidarf@mcmaster.ca',
    package_dir = {'cdft': 'cdft'},
    packages=['cdft', 'cdft.tools', 'cdft.tools.tests',
              'cdft.analyze', 'cdft.analyze.tests'],
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2',
        'Topic :: Science/Engineering :: Molecular Science'
        ],
     )
