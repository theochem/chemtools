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
"""The Input-Output (IO) Module that is user friendly."""


def make_molecule(*args, **kwargs):
    """Creates a Molecule instance using the given package.

    Parameters
    ----------
    args : iterable
        Arguments that will be used to create the Molecule instance
        If only one argument is given, it will be assumed that the argument is a filename and the
        classmethod from_file will be used instead of the initializer
    kwargs : dict
        Keyword arguments that will be used to create the Molecule instance
        package_name : str
            Name of the package that will be used to create the Molecule instance
            The wrapper for this package must exist in ChemTools
            By default,
        For others, see appropriate Molecule class for details.

    Returns
    -------
    molecule : Molecule
        Appropriate instance of Molecule

    Raises
    ------
    NotImplementedError
        If there are no packages available that can be used with ChemTools
        If the specified package cannot be found
    """
    if 'package_name' not in kwargs:
        for abs_classname in ['chemtools.utils.wrappers.Psi4Molecule',
                              'chemtools.utils.wrappers.HortonMolecule',]:
            try:
                modulename, classname = abs_classname.rsplit('.', 1)
                Molecule = getattr(__import__(modulename), classname)
                if len(args) == 1 and isinstance(args[0], str):
                    return Molecule.from_file(*args)
                else:
                    return Molecule(*args, **kwargs)
            except (AttributeError, ImportError):
                pass
        else:
            raise NotImplementedError('Cannot find any of the packages compatible with ChemTools.')

    package_name = kwargs['package_name'].lower()
    if package_name == 'horton':
        from chemtools.utils.wrappers import HortonMolecule
        return HortonMolecule(*args, **kwargs)
    else:
        raise NotImplementedError('Given package, {0}, is not supported with '
                                  'ChemTools.'.format(package_name))
