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

import os
from distutils.command.install_data import install_data
from distutils.core import setup
from glob import glob


class my_install_data(install_data):
    """Add a datadir.txt file that points to the root for the data files. It is
       otherwise impossible to figure out the location of these data files at
       runtime.
    """

    def run(self):
        # Do the normal install_data
        install_data.run(self)
        # Create the file datadir.txt. It's exact content is only known
        # at installation time. By default, it is the installation prefix
        # passed to setup.py, but one can override it using the env var
        # INSTALL_DATA, which may be useful for packaging, or any other
        # situation where the installed files are moved to a new location
        # afterwards.
        my_install_dir = os.getenv("INSTALL_DIR", self.install_dir)
        # Loop over all packages in this project and write the data_dir.txt
        # file only in the main package. Usualy, there is only one that matters.
        dist = self.distribution
        libdir = dist.command_obj["install_lib"].install_dir
        for name in dist.packages:
            # If a package contains a dot, e.g. horton.test, then don't write
            # the file data_dir.txt.
            if '.' not in name:
                destination = os.path.join(libdir, name, "data_dir.txt")
                print "Creating %s" % destination
                if not self.dry_run:
                    with open(destination, "w") as f:
                        print >> f, my_install_dir


setup(
    name='chemtools',
    version='0.9.0',
    description='Package of Chemical Tools for Interpreting Quantum Chemistry Calculations',
    author='Ayers Group',
    author_email='horton.chemtools@gmail.com',
    package_dir={'chemtools': 'chemtools'},
    packages=[
        'chemtools', 'chemtools.toolbox', 'chemtools.conceptual',
        'chemtools.orbtools', 'chemtools.denstools', 'chemtools.utils',
    ],
    scripts=glob("scripts/*.py"),
    cmdclass={
        'install_data': my_install_data,
    },
    data_files=[
        ('share/chemtools', glob('data/*.*')),
        ('share/chemtools/test', glob('data/test/*.*')),
        ('share/chemtools/examples', glob('data/examples/*.*')),
    ],
    classifiers=[
        'Environment :: Console', 'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2',
        'Topic :: Science/Engineering :: Molecular Science'
    ],
    requires=[
        'numpy', 'horton', 'numpy', 'sphinx', 'matplotlib', 'PIL', 'mayavi',
        'Image', 'sympy', 'scipy'
    ])
