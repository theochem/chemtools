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
"""
Module for Conceptual Density Functional Theory Analysis of Quantum Chemistry Output Files.

This modules contains wrappers which take outputs of quantum chemistry software and
compute various conceptual density functional theory (DFT) descriptive tools.
"""

import sys
import logging

from chemtools.toolbox.utils import check_arg_molecule, get_matching_attr
from chemtools.toolbox.utils import get_dict_energy, get_dict_density, get_dict_population
from chemtools.conceptual.linear import LinearGlobalTool, LinearLocalTool, LinearCondensedTool
from chemtools.conceptual.quadratic import QuadraticGlobalTool, QuadraticLocalTool
from chemtools.conceptual.quadratic import QuadraticCondensedTool
from chemtools.conceptual.exponential import ExponentialGlobalTool
from chemtools.conceptual.rational import RationalGlobalTool
from chemtools.conceptual.general import GeneralGlobalTool

if sys.version_info.major == 2:
    from chemtools.wrappers2.molecule import Molecule
    from chemtools.wrappers2.grid import MolecularGrid
else:
    from chemtools.wrappers3.molecule import Molecule
    from chemtools.wrappers3.grid import MolecularGrid

try:
    from pathlib2 import Path
except ImportError:
    from pathlib import Path


class BaseConceptualDFT(object):
    """Base class for conceptual density functional theory (DFT) analysis."""

    def __init__(self, dict_values, dict_models, model, coordinates, numbers, **kwargs):
        r"""
        Initialize class.

        Parameters
        ----------
        dict_values : dict
            Dictionary of number_electrons :math:`N` (keys) and corresponding property (values).
        dict_models : dict
            Dictionary of energy_models (keys) and corresponding model_object (values).
        model : str
            Energy model.
        coordinates : ndarray
            Coordinates of atomic centers.
        numbers : ndarray
            Atomic number of atomic centers.

        """
        # check model
        if model.lower() not in dict_models.keys():
            raise ValueError("Model={0} is not available!".format(model.lower()))
        self._model = model.lower()

        # check shape of coordinates
        if coordinates is not None and (coordinates.ndim != 2 or coordinates.shape[1] != 3):
            raise ValueError("Argument coordinate should be a 2D-array with 3 columns! "
                             "Given coordinates.shape={0}".format(coordinates.shape))
        # check number of atoms given by numbers and coordinates match
        if numbers is not None and coordinates is not None and len(numbers) != len(coordinates):
            raise ValueError("Numbers & coordinates should represent same number of atoms! "
                             "{0}!={1}".format(len(numbers), len(coordinates)))
        self._coordinates = coordinates
        self._numbers = numbers

        if self._model == "general":
            # self._tool = select_tool[model](*args, **kwargs)
            raise NotImplementedError("Model={0} is not covered yet!".format(self._model))

        # check number of electrons are integers
        if not all([isinstance(item, (int, float)) for item in dict_values.keys()]):
            raise ValueError("For model={0}, integer number of electrons are required! "
                             "#electrons={1}".format(self._model, dict_values.keys()))

        # make an instance of global tool
        self._tool = dict_models[self._model](dict_values, **kwargs)

        # print screen information
        self._log_init()

    def __getattr__(self, attr):
        """Return class attribute."""
        value = getattr(self._tool, attr, "error")
        if isinstance(value, str) and value == "error":
            raise AttributeError("Attribute {0} does not exist!".format(attr))
        return value

    @property
    def model(self):
        """Energy model."""
        return self._model

    @property
    def coordinates(self):
        """Cartesian coordinates of atoms."""
        return self._coordinates

    @property
    def numbers(self):
        """Atomic numbers of atoms."""
        return self._numbers

    def _log_init(self):
        """Print an initial informative message."""
        logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
        logging.info("Initialize   : {0}".format(self.__class__))
        logging.info("Energy Model : {0}".format(self.model))
        # logging.info("Reference N0  : {0}".format(self._tool.n0))

    @staticmethod
    def load_file(fnames):
        """Return `Molecule` instances corresponding to fnames.

        Parameters
        ----------
        fnames : str or Sequence of str
            Strings specifying the path to molecule's file, or sequence of strings specifying
            path to molecule files.
        """
        if isinstance(fnames, (str, unicode, Path)):
            # case of one file not given as a list
            molecule = Molecule.from_file(fnames)
        elif len(fnames) == 1 and isinstance(fnames[0], (str, unicode)):
            # case of one file given as a list
            molecule = Molecule.from_file(fnames[0])
        else:
            # case of multiple files
            molecule = [Molecule.from_file(fname) for fname in fnames]
        return molecule


class GlobalConceptualDFT(BaseConceptualDFT):
    r"""
    Global conceptual density functional theory (DFT) analysis of quantum chemistry output files.

    Global conceptual DFT reactivity descriptors assign a single property value to the entire
    molecule, :math:`P_{\text{global}}`, which describes the intrinsic reactivity of the molecule
    as a whole.

    Note: If only one molecule is provided, the frontier molecular orbital (FMO) approach is
          invoked, otherwise finite difference (FD) approach is taken.
    """

    def __init__(self, dict_values, model="quadratic", coordinates=None, numbers=None):
        r"""
        Initialize class.

        Parameters
        ----------
        dict_values : dict
            Dictionary of number of electrons :math:`N` (key) and corresponding energy values
            :math:`E(N)` (value).
        model : str, optional
            Energy model used to calculate global reactivity descriptors.
            The available models include: "linear", "quadratic", "exponential", "rational",
            and "general". Please see "" for more information.
        coordinates : ndarray, optional
            Atomic coordinates of atomic centers.
        numbers : ndarray, optional
            Atomic number of atomic centers.

        """
        # available models for global tools
        dict_models = {"linear": LinearGlobalTool, "quadratic": QuadraticGlobalTool,
                       "exponential": ExponentialGlobalTool, "rational": RationalGlobalTool,
                       "general": GeneralGlobalTool}
        super(GlobalConceptualDFT, self).__init__(dict_values, dict_models, model,
                                                  coordinates, numbers)

    def __repr__(self):
        """Print table of available class attributes and methods."""
        # get list of available attributes and methods
        avs = dir(self._tool)
        # get sorted list of public attributes
        attrs = [atr for atr in avs if not callable(getattr(self, atr)) and not atr.startswith("_")]
        attrs.remove("params")
        attrs.sort()
        # get sorted list of public methods
        methods = [atr for atr in avs if callable(getattr(self, atr)) and not atr.startswith("_")]
        methods.sort()
        content = "\nAttributes of {0} global model:\n{1}\n".format(self._model, "-" * 50)
        for attr in attrs:
            if attr == 'dict_energy':
                continue
            value = getattr(self._tool, attr)
            if value is not None and not hasattr(value, "__iter__"):
                content += "\n%s   % .6f" % (attr.ljust(25), value)
            elif value is not None and hasattr(value, "__iter__"):
                content += "\n%s   [%s]" % (attr.ljust(25), ", ".join(["% .6f" % v for v in value]))
            else:
                content += "\n%s   %s" % (attr.ljust(25), " ---")
        content += "\n\nMethods of {0} global model:\n{1}\n".format(self._model, "-" * 50)
        content += "\n".join(methods) + "\n"
        return content

    @classmethod
    def from_file(cls, fname, model):
        r"""
        Initialize class from calculation output file(s).

        Parameters
        ----------
        fname : str, sequence of str
            String specifying the path to a molecule's file, or sequence of strings specifying
            path to molecule files.
        model : str
            Energy model used to calculate global reactivity descriptors.
        """
        molecules = cls.load_file(fname)
        return cls.from_molecule(molecules, model)

    @classmethod
    def from_molecule(cls, molecule, model):
        r"""
        Initialize class from `Molecule` object(s).

        Parameters
        ----------
        molecule : Molecule or Sequence of Molecule
            Instance of Molecule class, or sequence of Molecule class instances.
        model : str
            Energy model used to calculate global reactivity descriptors.
        """
        # check molecule
        molecule = check_arg_molecule(molecule)
        # get atomic number and atomic coordinates
        number = get_matching_attr(molecule, "numbers", 1.e-8)
        coords = get_matching_attr(molecule, "coordinates", 1.e-4)
        dict_energy = get_dict_energy(molecule)
        return cls(dict_energy, model, coords, number)


class LocalConceptualDFT(BaseConceptualDFT):
    r"""
    Local conceptual density functional theory (DFT) analysis of quantum chemistry output files.

    Local conceptual DFT reactivity descriptors :math:`p_{\text{local}} (\mathbf{r})`, assign a
    value to every point :math:`\mathbf{r}` in space which provide insight into the reactivity of
    the molecule at point.

    Note: If only one molecule is provided, the frontier molecular orbital (FMO) approach is
          invoked, otherwise finite difference (FD) approach is taken. If FD approach is taken,
          the geometries and atomic numbers of all molecules need to be the same.
    """

    def __init__(self, dict_density, model="quadratic", coordinates=None, numbers=None,
                 n_max=None, global_softness=None):
        r"""
        Initialize class.

        Parameters
        ----------
        dict_density : dict
            Dictionary of number of electrons :math:`N` (key) and corresponding density array
            :math:`\rho_N(\mathbf{r})` (value).
        model : str, optional
            Energy model used to calculate local reactivity descriptors.
            The available models include: "linear", "quadratic", and "general".
            Please see "" for more information.
        coordinates : ndarray, optional
            Atomic coordinates of atomic centers.
        numbers : ndarray, optional
            Atomic number of atomic centers.
        n_max : float, optional
            Maximum number of electrons that system can accept, i.e. :math:`N_{\text{max}}`.
            See :attr:`BaseGlobalTool.n_max`.
        global_softness : float, optional
            Global softness. See :attr:`BaseGlobalTool.softness`.

        """
        # available models for local tools
        dict_models = {"linear": LinearLocalTool, "quadratic": QuadraticLocalTool}
        # check density array shape
        for index, value in enumerate(dict_density.values()):
            if index == 0:
                shape = value.shape
            elif value.shape != shape:
                raise ValueError("Argument dict_density should have density arrays (values) "
                                 "with the same shape! {0} != {1}".format(value.shape, shape))
        kwargs = {"n_max": n_max, "global_softness": global_softness}
        super(LocalConceptualDFT, self).__init__(dict_density, dict_models, model,
                                                 coordinates, numbers, **kwargs)

    def __repr__(self):
        """Print table of available class attributes and methods."""
        # get list of available attributes and methods
        avs = dir(self._tool) + dir(self)
        # get sorted list of public attributes
        attrs = [atr for atr in avs if not callable(getattr(self, atr)) and not atr.startswith("_")]
        attrs.sort()
        # remove n0, because it is both an attr of self._tool and self (duplicate)
        # attrs.remove("n0")
        # get sorted list of public methods
        methods = [atr for atr in avs if callable(getattr(self, atr)) and not atr.startswith("_")]
        methods.sort()
        content = "Attributes of {0} local model:\n{1}\n".format(self._model, "-" * 50)
        content += "\n".join(attrs)
        content += "\n\nMethods of {0} local model:\n{1}\n".format(self._model, "-" * 50)
        content += "\n".join(methods) + "\n"
        return content

    @classmethod
    def from_file(cls, fname, model, points):
        r"""
        Initialize class from calculation output file(s).

        Parameters
        ----------
        fname : str, sequence of str
            String specifying the path to a molecule's file, or sequence of strings specifying
            path to molecule files.
        model : str
            Energy model used to calculate local reactivity descriptors.
            Available models are "linear" and "quadratic".
        points : np.array
            Coordinates of points on which the local properties are evaluated given as a 2D
            array with 3 columns.
        """
        molecules = cls.load_file(fname)
        return cls.from_molecule(molecules, model, points)

    @classmethod
    def from_molecule(cls, molecule, model, points):
        r"""
        Initialize class from `Molecule` object(s).

        Parameters
        ----------
        molecule : Molecule or Sequence of Molecule
            Instance of Molecule class, or sequence of Molecule class instances.
        model : str
            Energy model used to calculate local reactivity descriptors.
            Available models are "linear" and "quadratic".
        points : np.array
            Coordinates of points on which the local properties are evaluated given as a 2D
            array with 3 columns.
        """
        # check molecule
        molecule = check_arg_molecule(molecule)
        # if points is None:
        #     # make molecular grid which also checks for matching atomic numbers & coordinates
        #     grid = get_molecular_grid(molecule, grid=None)
        #     points, numbers, coords = grid.points, grid.numbers, grid.centers
        # get atomic number, coordinates, N_max, and global softness
        gcdft = GlobalConceptualDFT.from_molecule(molecule, model)
        coords, nums, n_max, s = gcdft.coordinates, gcdft.numbers, gcdft.n_max, gcdft.softness
        dict_dens = get_dict_density(molecule, points)
        return cls(dict_dens, model, coords, nums, n_max, s)


class CondensedConceptualDFT(BaseConceptualDFT):
    r"""
    Condensed conceptual density functional theory (DFT) analysis of quantum chemistry output files.

    Condensed conceptual DFT reactivity descriptors :math:`\{p_A\}_{A=1}^{N_A}`, partition the
    local descriptors between atoms in molecules, and integrate it to obtain condense property.

    Note: If only one molecule is provided, the frontier molecular orbital (FMO) approach is
          invoked, otherwise finite difference (FD) approach is taken. If FD approach is taken,
          the geometries and atomic numbers of all molecules can be different only for fragment
          of molecular response (approach='RMF').
    """

    def __init__(self, dict_population, model="quadratic", coordinates=None, numbers=None,
                 n_max=None, global_softness=None):
        r"""
        Initialize class.

        Parameters
        ----------
        dict_population : dict
            Dictionary of number of electrons :math:`N` (key) and corresponding atomic population
            array (value).
        model : str, optional
            Energy model used to calculate condensed reactivity descriptors.
            The available models include: "linear", "quadratic", and "general".
            Please see "" for more information.
        coordinates : ndarray, optional
            Atomic coordinates of atomic centers.
        numbers : ndarray, optional
            Atomic number of atomic centers.
        n_max : float, optional
            Maximum number of electrons that system can accept, i.e. :math:`N_{\text{max}}`.
            See :attr:`BaseGlobalTool.n_max`.
        global_softness : float, optional
            Global softness. See :attr:`BaseGlobalTool.softness`.

        """
        # available models for local tools
        dict_models = {"linear": LinearCondensedTool, "quadratic": QuadraticCondensedTool}
        kwargs = {"n_max": n_max, "global_softness": global_softness}
        super(CondensedConceptualDFT, self).__init__(dict_population, dict_models, model,
                                                     coordinates, numbers, **kwargs)

    def __repr__(self):
        """Print table of available class attributes and methods."""
        # get list of available attributes and methods
        avs = dir(self._tool)
        # get sorted list of public attributes
        attrs = [atr for atr in avs if not callable(getattr(self, atr)) and not atr.startswith("_")]
        attrs.sort()
        attrs.remove("n_ref")
        # get sorted list of public methods
        methods = [atr for atr in avs if callable(getattr(self, atr)) and not atr.startswith("_")]
        content = "Attributes of {0} condensed model:\n{1}\n".format(self._model, "-" * 50)
        content += "\n".join(attrs)
        content += "\n\nMethods in {0} condensed model:\n{1}\n".format(self._model, "-" * 50)
        content += "\n".join(methods) + "\n"
        return content

    @classmethod
    def from_file(cls, fname, model, approach="FMR", scheme="h", **kwargs):
        r"""
        Initialize class from calculation output file(s).

        Parameters
        ----------
        fname : str, sequence of str
            String specifying the path to a molecule's file, or sequence of strings specifying
            path to molecule files.
        model : str
            Energy model used to calculate condensed reactivity descriptors.
            Available models are "linear" and "quadratic".
        approach : str, optional
            Choose between "FMR" (fragment of molecular response) or "RMF"
            (response of molecular fragment).
        scheme : str, optional
            Partitioning scheme. Options: "h", "hi", "mbis".
        kwargs : dict, optional
            Extra keyword arguments required for partitioning, like 'grid' and 'proatomdb'.
        """
        molecules = cls.load_file(fname)
        return cls.from_molecule(molecules, model, approach, scheme, **kwargs)

    @classmethod
    def from_molecule(cls, molecule, model, approach="FMR", scheme="h", **kwargs):
        r"""
        Initialize class from `Molecule` object(s).

        Parameters
        ----------
        molecule : Molecule or Sequence of Molecule
            Instance of Molecule class, or sequence of Molecule class instances.
        model : str
            Energy model used to calculate condensed reactivity descriptors.
            Available models are "linear" and "quadratic".
        approach : str, optional
            Choose between "FMR" (fragment of molecular response) or "RMF"
            (response of molecular fragment).
        scheme : str, optional
            Partitioning scheme. Options: "h", "hi", "mbis".
        kwargs : dict, optional
            Extra keyword arguments required for partitioning, like 'grid' and 'proatomdb'.
        """
        # check molecule
        molecule = check_arg_molecule(molecule)
        # check type of grid
        if "grid" in kwargs.keys() and not isinstance(kwargs["grid"], MolecularGrid):
            raise ValueError("Currently, only 'MolecularGrid' is supported for condensing!")
        # get atomic number, coordinates, N_max, and global softness
        gcdft = GlobalConceptualDFT.from_molecule(molecule, model)
        coords, nums, n_max, s = gcdft.coordinates, gcdft.numbers, gcdft.n_max, gcdft.softness
        # get dictionary of populations
        dict_pops = get_dict_population(molecule, approach, scheme, **kwargs)
        return cls(dict_pops, model, coords, nums, n_max, s)
