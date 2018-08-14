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
"""
Module for Conceptual Density Functional Theory Analysis of Quantum Chemistry Output Files.

This modules contains wrappers which take outputs of quantum chemistry software and
compute various conceptual density functional theory (DFT) descriptive tools.
"""

import numpy as np

from horton import log
from horton import BeckeMolGrid, ProAtomDB
from horton.scripts.wpart import wpart_schemes

from chemtools.utils.molecule import BaseMolecule
from chemtools.toolbox.molecule import make_molecule
from chemtools.conceptual.linear import LinearGlobalTool, LinearLocalTool
from chemtools.conceptual.quadratic import QuadraticGlobalTool, QuadraticLocalTool
from chemtools.conceptual.exponential import ExponentialGlobalTool
from chemtools.conceptual.rational import RationalGlobalTool
from chemtools.conceptual.general import GeneralGlobalTool


__all__ = ["GlobalConceptualDFT", "LocalConceptualDFT", "CondensedConceptualDFT"]


class BaseConceptualDFT(object):
    """Base class for conceptual density functional theory (DFT) analysis."""

    def __init__(self, dict_values, dict_models, model, coordinates, numbers):
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
        self._tool = dict_models[self._model](dict_values)

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

    @staticmethod
    def load_file(file_names):
        """Return `BaseMolecule` instances corresponding to file_names.

        Parameters
        ----------
        file_names : str or Sequence of str
            Strings specifying the path to molecule's file, or sequence of strings specifying
            path to molecule files.
        """
        if isinstance(file_names, (str, unicode)):
            # case of one file not given as a list
            molecule = make_molecule(file_names)
        elif len(file_names) == 1 and isinstance(file_names[0], (str, unicode)):
            # case of one file given as a list
            molecule = make_molecule(file_names[0])
        else:
            # case of multiple files
            molecule = [make_molecule(filename) for filename in file_names]
        return molecule

    @staticmethod
    def _get_homo_lumo_data(molecule):
        r"""Return the HOMO and LUMO energy and HOMO and LUMO spin, respectively.

        Parameters
        ----------
        molecule : BaseMolecule
            Instance of BaseMolecule class.
        """
        # get homo & lumo energy and spin
        homo_energy = molecule.homo_energy[0]
        homo_spin = "a"
        if molecule.homo_energy[1] > homo_energy:
            homo_energy = molecule.homo_energy[1]
            homo_spin = "b"
        lumo_energy = molecule.lumo_energy[0]
        lumo_spin = "a"
        if molecule.lumo_energy[1] < lumo_energy:
            lumo_energy = molecule.lumo_energy[1]
            lumo_spin = "b"
        return homo_energy, lumo_energy, homo_spin, lumo_spin

    @staticmethod
    def _get_matching_attr(molecule, attr, accuracy=1.e-6):
        r"""Return the specified attribute after checking that it matches between molecules.

        Parameters
        ----------
        molecule : BaseMolecule or Sequence of BaseMolecule
            Instance of BaseMolecule class, or sequence of BaseMolecule class instances.
        """
        if isinstance(molecule, BaseMolecule):
            # get attribute for single molecule
            ref = getattr(molecule, attr)
        elif np.all([isinstance(mol, BaseMolecule) for mol in molecule]):
            # check whether attr matches between molecules
            for index, mol in enumerate(molecule):
                if index == 0:
                    ref = getattr(mol, attr)
                    continue
                item = getattr(mol, attr)
                if item.shape != ref.shape or not np.max(abs(ref - item)) < accuracy:
                    raise ValueError(
                        "Molecule 1 & {0} have different {1}!".format(index + 1, attr))
        else:
            raise ValueError("Argument molecule not recognized!")
        return ref

    @staticmethod
    def get_dict_energy(molecule):
        r"""Return dictionary of number of electrons and corresponding energy values.

        Parameters
        ----------
        molecule : BaseMolecule or Sequence of BaseMolecule
            Instance of BaseMolecule class, or sequence of BaseMolecule class instances.
            In the case of one molecule, the Frontier Orbital Molecule (FMO) approach is used
            to get the energy values of :math:`E(N + 1)` and :math:`E(N - 1)`.
        """
        if isinstance(molecule, BaseMolecule):
            # get homo/lumo energy
            homo_e, lumo_e, _, _ = BaseConceptualDFT._get_homo_lumo_data(molecule)
            nelec = sum(molecule.nelectrons)
            # store number of electron and energy in a dictionary
            energies = {nelec: molecule.energy,
                        nelec + 1: molecule.energy + lumo_e,
                        nelec - 1: molecule.energy - homo_e}
        elif np.all([isinstance(mol, BaseMolecule) for mol in molecule]):
            # store number of electron and energy in a dictionary
            energies = {}
            for mol in molecule:
                # get number of electrons & check that it doesn't already exist
                nelec = sum(mol.nelectrons)
                if nelec in energies.keys():
                    raise ValueError("Two molecules have {0} electrons!".format(nelec))
                # store number of electrons and energy in a dictionary
                energies[nelec] = mol.energy
        else:
            raise ValueError("Argument molecule not recognized!")
        return energies

    @staticmethod
    def get_dict_density(molecule, points):
        r"""Return dictionary of number of electrons and corresponding density values.

        Parameters
        ----------
        molecule : BaseMolecule or Sequence of BaseMolecule
            Instance of BaseMolecule class, or sequence of BaseMolecule class instances.
            In the case of one molecule, the Frontier Orbital Molecule (FMO) approach is used
            to get the density of :math:`\rho_{N + 1}(\mathbf{r})` and
            :math:`\rho_{N - 1}(\mathbf{r})`.
        points : ndarray
           The 2D array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        """
        if isinstance(molecule, BaseMolecule):
            # get homo/lumo energy and spin
            _, _, homo_s, lumo_s = BaseConceptualDFT._get_homo_lumo_data(molecule)
            # compute homo & lumo density
            spin_to_index = {"a": 0, "b": 1}
            homo_dens = molecule.compute_density(points, homo_s,
                                                 molecule.homo_index[spin_to_index[homo_s]])
            lumo_dens = molecule.compute_density(points, lumo_s,
                                                 molecule.lumo_index[spin_to_index[lumo_s]])
            # store number of electron and density in a dictionary
            nelec = sum(molecule.nelectrons)
            dens = molecule.compute_density(points, "ab", None)
            densities = {nelec: dens,
                         nelec + 1: dens + lumo_dens,
                         nelec - 1: dens - homo_dens}
        elif np.all([isinstance(mol, BaseMolecule) for mol in molecule]):
            # compute and record densities on given points in a dictionary
            densities = {}
            for mol in molecule:
                # get number of electrons
                nelec = sum(mol.nelectrons)
                if nelec in densities.keys():
                    raise ValueError("Two molecules have {0} electrons!".format(nelec))
                # store density for a given number of electron
                densities[nelec] = mol.compute_density(points, "ab", None)
        else:
            raise ValueError("Argument molecule not recognized!")
        return densities

    @staticmethod
    def get_dict_population(molecule, approach, grid, wpart, kwargs):
        r"""Return dictionary of number of electrons and corresponding atomic charges values.

        Parameters
        ----------
        molecule : BaseMolecule or Sequence of BaseMolecule
            Instance of BaseMolecule class, or sequence of BaseMolecule class instances.
        approach : str, optional
            Choose between "FMR" (fragment of molecular response) or "RMF"
            (response of molecular fragment).
        grid: BeckeMolGrid, optional
            Molecular grid used for partitioning.
        wpart :
        kwargs :
        """
        # check whether molecules have the same coordinates
        try:
            BaseConceptualDFT._get_matching_attr(molecule, "coordinates", 1.e-4)
            same_coordinates = True
        except ValueError:
            if approach.lower() == "fmr":
                raise ValueError("When geometries of molecules are different, only approach='RMF' "
                                 "is possible! Given approach={0}".format(approach.upper()))
            same_coordinates = False

        # check or generate molecular grid
        grid = BaseConceptualDFT._get_molecular_grid(molecule, grid)
        # compute dictionary of number of electron and density
        dict_dens = BaseConceptualDFT.get_dict_density(molecule, grid.points)

        # find reference molecule
        if isinstance(molecule, BaseMolecule):
            mol0 = molecule
        elif len(molecule) == 1 and isinstance(molecule[0], BaseMolecule):
            mol0 = molecule[0]
        elif np.all([isinstance(mol, BaseMolecule) for mol in molecule]):
            if len(molecule) != 3:
                raise ValueError("Conceptual DFT within FD approach, currently works for "
                                 "only 3 molecules! Given {0} molecules.".format(len(molecule)))

            # reference molecule is the middle molecule (for 3 molecules)
            dict_mols = {sum(mol.nelectrons): mol for mol in molecule}
            mol0 = dict_mols.pop(sorted(dict_mols.keys())[1])

        # compute population of reference molecule
        part0 = wpart(mol0.coordinates, mol0.numbers, mol0.pseudo_numbers, grid,
                      dict_dens[sum(mol0.nelectrons)], **kwargs)
        part0.do_all()
        # record population of reference system
        dict_pops = dict([(sum(mol0.nelectrons), part0["populations"])])
        del dict_dens[sum(mol0.nelectrons)]

        # compute and record populations given grid in a dictionary
        for nelec, dens in dict_dens.iteritems():
            # make sure there is no repetition
            if nelec in dict_pops.keys():
                raise ValueError("")
            if approach.lower() == "fmr":
                # fragment of molecular response
                pops = CondensedConceptualDFT.condense_to_atoms(dens, part0)
            elif approach.lower() == "rmf":
                # response of molecular fragment
                if not same_coordinates:
                    mol = dict_mols[nelec]
                    grid = BaseConceptualDFT._get_molecular_grid(molecule, None)
                    dens = molecule.compute_density(grid.points, "ab", None)
                parts = wpart(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                              grid, dens, **kwargs)

                parts.do_all()
                pops = parts["populations"]
            else:
                raise ValueError("Partitioning approach {0} is not recognized!".format(approach))
            # Store number of electron and populations in a dictionary
            dict_pops[nelec] = pops

        return dict_pops

    @staticmethod
    def _get_molecular_grid(molecule, grid=None):
        r"""Return molecular grid or check that given grid is consistent with molecule.

        Parameters
        ----------
        molecule : BaseMolecule or Sequence of BaseMolecule
            Instance of BaseMolecule class, or sequence of BaseMolecule class instances.
        grid : BeckeMolGrid, optional
            Instance or BeckeMolGrid. If `None`, a default `BeckeMolGrid` is returned.
        """
        # check type of grid
        if grid is not None and not isinstance(grid, BeckeMolGrid):
            raise ValueError("Currently, only 'BeckeMolGrid' is supported for condensing!")
        # check grid or make grid
        if grid is not None and isinstance(molecule, BaseMolecule):
            # check coordinates and atomic numbers of grid and molecule match
            if not np.all(abs(grid.centers - molecule.coordinates) < 1.e-4):
                raise ValueError("Coordinates of grid and molecule do not match!")
            if not np.all(abs(grid.numbers - molecule.numbers) < 1.e-4):
                raise ValueError("Atomic numbers of grid and molecule do not match!")
        elif grid is not None and all([isinstance(mol, BaseMolecule) for mol in molecule]):
            for index, mol in enumerate(molecule):
                # check atomic numbers of grid and molecules match
                ref, numbers = grid.numbers, mol.numbers
                if ref.shape != numbers.shape or not np.max(abs(ref - numbers)) < 1.e-6:
                    raise ValueError("Atomic numbers of grid and molecule {0} "
                                     "do not match!".format(index + 1))
                # check coordinates of grid and molecules match
                ref, coord = grid.centers, mol.coordinates
                if ref.shape != coord.shape or not np.max(abs(ref - coord)) < 1.e-4:
                    raise ValueError("Coordinates of grid and molecule {0} "
                                     "do not match!".format(index + 1))
        else:
            # make default grid
            number = BaseConceptualDFT._get_matching_attr(molecule, "numbers", 1.e-8)
            pseudo = BaseConceptualDFT._get_matching_attr(molecule, "pseudo_numbers", 1.e-8)
            coords = BaseConceptualDFT._get_matching_attr(molecule, "coordinates", 1.e-4)
            grid = BeckeMolGrid(coords, number, pseudo, agspec="fine", random_rotate=False,
                                mode="keep")
        return grid

    @staticmethod
    def _get_part_specifications(scheme, proatomdb, numbers, kwargs):
        """
        """
        # get partitioning class
        if scheme.lower() not in wpart_schemes:
            raise ValueError("Partitioning scheme={0} not supported! "
                             "Select from: {1}".format(scheme, wpart_schemes.keys()))
        wpart = wpart_schemes[scheme]
        # make proatom database
        # kwargs = {}
        if scheme.lower() not in ["mbis", "b"]:
            if proatomdb is None:
                proatomdb = ProAtomDB.from_refatoms(numbers)
            kwargs["proatomdb"] = proatomdb
        return wpart, kwargs


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
        content = "\nAvailable attributes in {0} global model:\n{1}\n".format(self._model, "-" * 50)
        for attr in attrs:
            value = getattr(self._tool, attr)
            if value is not None and not hasattr(value, "__iter__"):
                content += "\n%s   % .6f" % (attr.ljust(25), value)
            elif value is not None and hasattr(value, "__iter__"):
                content += "\n%s   [%s]" % (attr.ljust(25), ", ".join(["% .6f" % v for v in value]))
            else:
                content += "\n%s   %s" % (attr.ljust(25), " ---")
        content += "\n\nAvailable methods in {0} global model:\n{1}\n".format(self._model, "-" * 50)
        content += "\n".join(methods) + "\n"
        return content

    def _log_init(self):
        """Print an initial informative message."""
        if log.do_medium:
            log("Initialize: %s" % self.__class__)
            log.deflist([("Energy Model", self._model),
                         ("Parameters", self._tool.params),
                         ("Reference Energy", self._tool.dict_energy)])
            log.blank()

    @classmethod
    def from_file(cls, file_name, model):
        r"""
        Initialize class from calculation output file(s).

        Parameters
        ----------
        file_name : str, sequence of str
            String specifying the path to a molecule's file, or sequence of strings specifying
            path to molecule files.
        model : str
            Energy model used to calculate global reactivity descriptors.
        """
        molecules = cls.load_file(file_name)
        return cls.from_molecule(molecules, model)

    @classmethod
    def from_molecule(cls, molecule, model):
        r"""
        Initialize class from `BaseMolecule` object(s).

        Parameters
        ----------
        molecule : BaseMolecule or Sequence of BaseMolecule
            Instance of BaseMolecule class, or sequence of BaseMolecule class instances.
        model : str
            Energy model used to calculate global reactivity descriptors.
        """
        # get atomic number and atomic coordinates
        number = cls._get_matching_attr(molecule, "numbers", 1.e-8)
        coords = cls._get_matching_attr(molecule, "coordinates", 1.e-4)
        dict_energy = cls.get_dict_energy(molecule)
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

    def __init__(self, dict_density, model="quadratic", coordinates=None, numbers=None):
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
        super(LocalConceptualDFT, self).__init__(dict_density, dict_models, model,
                                                 coordinates, numbers)

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
        content = "Available attributes in {0} global model:\n{1}\n".format(self._model, "-" * 50)
        content += "\n".join(attrs)
        content += "\n\nAvailable methods in {0} global model:\n{1}\n".format(self._model, "-" * 50)
        content += "\n".join(methods) + "\n"
        return content

    def _log_init(self):
        """Print an initial informative message."""
        if log.do_medium:
            log("Initialize: %s" % self.__class__)
            log.deflist([("Energy Model", self._model),
                         ("Reference #Electrons", self._tool.n0)])
            log.blank()

    @classmethod
    def from_file(cls, file_name, model, points=None):
        r"""
        Initialize class from calculation output file(s).

        Parameters
        ----------
        file_name : str, sequence of str
            String specifying the path to a molecule's file, or sequence of strings specifying
            path to molecule files.
        model : str
            Energy model used to calculate local reactivity descriptors.
            Available models are "linear" and "quadratic".
        points : np.array, optional
            Points on which the local properties are evaluated. If `None`, a default BeckeMolGrid
            is generated to get the points.
        """
        molecules = cls.load_file(file_name)
        return cls.from_molecule(molecules, model, points)

    @classmethod
    def from_molecule(cls, molecule, model, points=None):
        r"""
        Initialize class from `BaseMolecule` object(s).

        Parameters
        ----------
        molecule : BaseMolecule or Sequence of BaseMolecule
            Instance of BaseMolecule class, or sequence of BaseMolecule class instances.
        model : str
            Energy model used to calculate local reactivity descriptors.
            Available models are "linear" and "quadratic".
        points : np.array, optional
            Points on which the local properties are evaluated. If `None`, a default BeckeMolGrid
            is generated to get the points.
        """
        if points is None:
            # make molecular grid which also checks for matching atomic numbers & coordinates
            grid = cls._get_molecular_grid(molecule, grid=None)
            points, numbers, coords = grid.points, grid.numbers, grid.centers
        else:
            numbers = cls._get_matching_attr(molecule, "numbers", 1.e-8)
            coords = cls._get_matching_attr(molecule, "coordinates", 1.e-4)
        dict_dens = cls.get_dict_density(molecule, points)
        return cls(dict_dens, model, coords, numbers)


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

    def __init__(self, dict_population, model="quadratic", coordinates=None, numbers=None):
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

        """
        # available models for local tools
        dict_models = {"linear": LinearLocalTool, "quadratic": QuadraticLocalTool}
        super(CondensedConceptualDFT, self).__init__(dict_population, dict_models, model,
                                                     coordinates, numbers)

    def __repr__(self):
        """Print table of available class attributes and methods."""
        # get list of available attributes and methods
        avs = dir(self._tool)
        # get sorted list of public attributes
        attrs = [atr for atr in avs if not callable(getattr(self, atr)) and not atr.startswith("_")]
        attrs.sort()
        attrs.remove("n0")
        # get sorted list of public methods
        methods = [atr for atr in avs if callable(getattr(self, atr)) and not atr.startswith("_")]
        content = "Available attributes in {0} global model:\n{1}\n".format(self._model, "-" * 50)
        content += "\n".join(attrs)
        content += "\n\nAvailable methods in {0} global model:\n{1}\n".format(self._model, "-" * 50)
        content += "\n".join(methods) + "\n"
        return content

    def _log_init(self):
        """Print an initial informative message."""
        if log.do_medium:
            log("Initialize: Condensed Class")
            # log("Initialize: %s" % self.__class__)
            log.deflist([("Energy Model", self._model),
                         ("Reference #Electrons", self._tool.n0)])
            log.blank()

    @classmethod
    def from_file(cls, file_name, model, approach="FMR", scheme="h", grid=None, proatomdb=None):
        r"""
        Initialize class from calculation output file(s).

        Parameters
        ----------
        file_name : str, sequence of str
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
        grid : BeckeMolGrid, optional
            Molecular grid used for partitioning.
        proatomdb: ProAtomDB
            Pro-atom database used for partitioning. Only "h" and "hi" requires that.
        """
        molecules = cls.load_file(file_name)
        return cls.from_molecule(molecules, model, approach, scheme, grid, proatomdb)

    @classmethod
    def from_molecule(cls, molecule, model, approach="FMR", scheme="h", grid=None,
                      proatomdb=None, **kwargs):
        r"""
        Initialize class from `BaseMolecule` object(s).

        Parameters
        ----------
        molecule : BaseMolecule or Sequence of BaseMolecule
            Instance of BaseMolecule class, or sequence of BaseMolecule class instances.
        model : str
            Energy model used to calculate condensed reactivity descriptors.
            Available models are "linear" and "quadratic".
        approach : str, optional
            Choose between "FMR" (fragment of molecular response) or "RMF"
            (response of molecular fragment).
        scheme : str, optional
            Partitioning scheme. Options: "h", "hi", "mbis".
        grid: BeckeMolGrid, optional
            Molecular grid used for partitioning.
        proatomdb : ProAtomDB
            Pro-atom database used for partitioning. Only "h" and "hi" requires that.
        kwargs :
        """
        numbers = cls._get_matching_attr(molecule, "numbers", 1.e-8)
        coords = cls._get_matching_attr(molecule, "coordinates", 1.e-4)
        # check and get partitioning object
        wpart, kwargs = cls._get_part_specifications(scheme, proatomdb, numbers, kwargs)
        # get dictionary of populations
        dict_pops = cls.get_dict_population(molecule, approach, grid, wpart, kwargs)
        return cls(dict_pops, model, coords, numbers)

    @staticmethod
    def condense_to_atoms(local_property, part):
        r"""
        Return condensed values of the local descriptor partitioned and integrated over atoms.

        Condense the local descriptor :math:`p_{\text{local}}\left(\mathbf{r}\right)` into
        atomic contributions :math:`\{P_A\}_{A=1}^N_{\text{atoms}}` defined as,

        .. math::
           P_A = \int \omega_A\left(\mathbf{r}\right)
                      p_{\text{local}}\left(\mathbf{r}\right) d\mathbf{r}

        Parameters
        ----------
        local_property : ndarray
            Local descriptor evaluated on grid.
        part : part instance
            Instance of `HORTON` partitioning calss.
        """
        condensed = np.zeros(part.natom)
        for index in range(part.natom):
            at_grid = part.get_grid(index)
            at_weight = part.cache.load("at_weights", index)
            wcor = part.get_wcor(index)
            local_prop = part.to_atomic_grid(index, local_property)
            condensed[index] = at_grid.integrate(at_weight, local_prop, wcor)
        return condensed
