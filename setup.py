#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2019 The ChemTools Development Team
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
# pragma pylint: disable=superfluous-parens

from distutils.core import setup

setup(
    name='chemtools',
    version='0.0.0',
    description='Package of Chemical Tools for Interpreting Quantum Chemistry Calculations',
    author='ChemTools Dev Team',
    author_email='horton.chemtools@gmail.com',
    package_dir={'chemtools': 'chemtools'},
    packages=[
        'chemtools',
        'chemtools.toolbox',
        'chemtools.conceptual',
        'chemtools.denstools',
        'chemtools.utils',
        'chemtools.wrappers',
        'chemtools.outputs',
        'chemtools.data',
        'chemtools.scripts',
        'chemtools.toolbox.test',
        'chemtools.conceptual.test',
        'chemtools.denstools.test',
        'chemtools.utils.test',
        'chemtools.wrappers.test',
        'chemtools.outputs.test',
    ],
    package_data={'chemtools.data': ['*.fchk', '*.cube', '*.wfn', '*.npz'],
                  'chemtools.data.examples': ['*.fchk', '*.cube', '*.wfn', '*.npz']},
    entry_points={
        'console_scripts': ['chemtools = chemtools.scripts.main:main'],
    },
    classifiers=[
        'Environment :: Console', 'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2',
        'Topic :: Science/Engineering :: Molecular Science'
    ],
    install_requires=[
        'numpy', 'matplotlib', 'Pillow', 'Image', 'sympy',
        'scipy', 'importlib_resources; python_version < "3.7"',
    ])
