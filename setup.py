#!/usr/bin/env python

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
