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
# pragma pylint: disable=too-many-branches,too-many-statements
"""Utility Functions of Toolbox Module."""


import sys
import numpy as np


if sys.version_info.major == 2:
    from chemtools.wrappers2.molecule import Molecule
    from chemtools.wrappers2.grid import MolecularGrid
    from horton import ProAtomDB
    from horton.scripts.wpart import wpart_schemes
else:
    from chemtools.wrappers3.molecule import Molecule
    from chemtools.wrappers3.grid import MolecularGrid


def check_arg_molecule(molecule):
    """Return molecule argument after checking.

    Parameters
    ----------
    molecule : Molecule or Sequence of Molecule
        Instance of Molecule class, or sequence of Molecule class instances.

    Returns
    -------
    molecule : Molecule or Sequence of Molecule
        Instance of Molecule or Sequence of Molecule with more than one instance.
    """
    if isinstance(molecule, Molecule):
        return molecule
    if hasattr(molecule, "__iter__") and len(molecule) == 1:
        # sequence of just one molecule
        return molecule[0]
    return molecule


def get_homo_lumo_data(molecule):
    r"""Return the HOMO and LUMO energy and HOMO and LUMO spin, respectively.

    Parameters
    ----------
    molecule : Molecule
        Instance of Molecule class.
    """
    # get homo & lumo energy and spin
    homo_energy = molecule.mo.homo_energy[0]
    homo_spin = "a"
    if molecule.mo.homo_energy[1] > homo_energy:
        homo_energy = molecule.mo.homo_energy[1]
        homo_spin = "b"
    lumo_energy = molecule.mo.lumo_energy[0]
    lumo_spin = "a"
    if molecule.mo.lumo_energy[1] < lumo_energy:
        lumo_energy = molecule.mo.lumo_energy[1]
        lumo_spin = "b"
    return homo_energy, lumo_energy, homo_spin, lumo_spin


def get_matching_attr(molecule, attr, accuracy=1.e-6):
    r"""Return the specified attribute after checking that it matches between molecules.

    Parameters
    ----------
    molecule : Molecule or Sequence of Molecule
        Instance of Molecule class, or sequence of Molecule class instances.
    attr : str
        The name ot attribute.
    accuracy : float, optional
        The accuracy for matching the attribute between different molecules.
    """
    if isinstance(molecule, Molecule):
        # get attribute for single molecule
        ref = getattr(molecule, attr)
    elif np.all([isinstance(mol, Molecule) for mol in molecule]):
        # check whether attr matches between molecules
        for index, mol in enumerate(molecule):
            if index == 0:
                ref = getattr(mol, attr)
                continue
            item = getattr(mol, attr)
            if item.shape != ref.shape or not np.max(abs(ref - item)) < accuracy:
                raise ValueError("Molecule 0 & {0} have different {1}!".format(index, attr))
    else:
        raise ValueError("Argument molecule not recognized!")
    return ref


def get_molecular_grid(molecule, grid=None):
    r"""Return molecular grid or check that given grid is consistent with molecule.

    Parameters
    ----------
    molecule : Molecule or Sequence of Molecule
        Instance of Molecule class, or sequence of Molecule class instances.
    grid : MolecularGrid, optional
        Instance or MolecularGrid. If `None`, a default `MolecularGrid` is returned.
    """
    # check grid or make grid
    if grid is not None and isinstance(molecule, Molecule):
        # check atomic numbers & coordinates of grid and molecule match
        ref, numbers = grid.numbers, molecule.numbers
        if ref.shape != numbers.shape or not np.max(abs(ref - numbers)) < 1.e-6:
            raise ValueError("Atomic numbers of grid and molecule do not match!")
        ref, coord = grid.centers, molecule.coordinates
        if ref.shape != coord.shape or not np.max(abs(ref - coord)) < 1.e-4:
            raise ValueError("Coordinates of grid and molecule do not match!")
    elif grid is not None and all([isinstance(mol, Molecule) for mol in molecule]):
        for index, mol in enumerate(molecule):
            # check atomic numbers of grid and molecules match
            ref, numbers = grid.numbers, mol.numbers
            if ref.shape != numbers.shape or not np.max(abs(ref - numbers)) < 1.e-6:
                raise ValueError("Atomic number of grid & molecule {0} do not match!".format(index))
            # check coordinates of grid and molecules match
            ref, coord = grid.centers, mol.coordinates
            if ref.shape != coord.shape or not np.max(abs(ref - coord)) < 1.e-4:
                raise ValueError("Coordinates of grid & molecule {0} do not match!".format(index))
    else:
        # make default grid
        number = get_matching_attr(molecule, "numbers", 1.e-8)
        pseudo = get_matching_attr(molecule, "pseudo_numbers", 1.e-8)
        coords = get_matching_attr(molecule, "coordinates", 1.e-4)
        grid = MolecularGrid(coords, number, pseudo, specs="insane", rotate=False)
    return grid


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


def get_dict_energy(molecule):
    r"""Return dictionary of number of electrons and corresponding energy values.

    Parameters
    ----------
    molecule : Molecule or Sequence of Molecule
        Instance of Molecule class, or sequence of Molecule class instances.
        In the case of one molecule, the Frontier Orbital Molecule (FMO) approach is used
        to get the energy values of :math:`E(N + 1)` and :math:`E(N - 1)`.
    """
    if isinstance(molecule, Molecule):
        # get homo/lumo energy
        homo_e, lumo_e, _, _ = get_homo_lumo_data(molecule)
        nelec = sum(molecule.mo.nelectrons)
        # store number of electron and energy in a dictionary
        energies = {nelec: molecule.energy,
                    nelec + 1: molecule.energy + lumo_e,
                    nelec - 1: molecule.energy - homo_e}
    elif np.all([isinstance(mol, Molecule) for mol in molecule]):
        # store number of electron and energy in a dictionary
        energies = {}
        for mol in molecule:
            # get number of electrons & check that it doesn't already exist
            nelec = sum(mol.mo.nelectrons)
            if nelec in energies.keys():
                raise ValueError("Two molecules have {0} electrons!".format(nelec))
            # store number of electrons and energy in a dictionary
            energies[nelec] = mol.energy
    else:
        raise ValueError("Argument molecule not recognized!")
    return energies


def get_dict_density(molecule, points):
    r"""Return dictionary of number of electrons and corresponding density values.

    Parameters
    ----------
    molecule : Molecule or Sequence of Molecule
        Instance of Molecule class, or sequence of Molecule class instances.
        In the case of one molecule, the Frontier Orbital Molecule (FMO) approach is used
        to get the density of :math:`\rho_{N + 1}(\mathbf{r})` and
        :math:`\rho_{N - 1}(\mathbf{r})`.
    points : ndarray
       The 2D array containing the cartesian coordinates of points on which density is
       evaluated. It has a shape (n, 3) where n is the number of points.
    """
    if isinstance(molecule, Molecule):
        # get homo/lumo energy and spin
        _, _, homo_s, lumo_s = get_homo_lumo_data(molecule)
        # compute homo & lumo density
        spin_to_index = {"a": 0, "b": 1}
        homo_dens = molecule.compute_density(points, homo_s,
                                             molecule.mo.homo_index[spin_to_index[homo_s]])
        lumo_dens = molecule.compute_density(points, lumo_s,
                                             molecule.mo.lumo_index[spin_to_index[lumo_s]])
        # store number of electron and density in a dictionary
        nelec = sum(molecule.mo.nelectrons)
        dens = molecule.compute_density(points, "ab", None)
        densities = {nelec: dens,
                     nelec + 1: dens + lumo_dens,
                     nelec - 1: dens - homo_dens}
    elif np.all([isinstance(mol, Molecule) for mol in molecule]):
        # compute and record densities on given points in a dictionary
        densities = {}
        for mol in molecule:
            # get number of electrons
            nelec = sum(mol.mo.nelectrons)
            if nelec in densities.keys():
                raise ValueError("Two molecules have {0} electrons!".format(nelec))
            # store density for a given number of electron
            densities[nelec] = mol.compute_density(points, "ab", None)
    else:
        raise ValueError("Argument molecule not recognized!")
    return densities


def get_dict_population(molecule, approach, scheme, **kwargs):
    r"""Return dictionary of number of electrons and corresponding atomic charges values.

    Parameters
    ----------
    molecule : Molecule or Sequence of Molecule
        Instance of Molecule class, or sequence of Molecule class instances.
    approach : str, optional
        Choose between "FMR" (fragment of molecular response) or "RMF"
        (response of molecular fragment).
    scheme : str
        Partitioning scheme.
    kwargs : optional
    """
    # check approach
    if approach.lower() not in ["rmf", "fmr"]:
        raise ValueError("Argument approach={0} is not valid.".format(approach))

    # case of populations available in molecule
    if scheme.lower() not in wpart_schemes:
        # check approach & molecule instances
        if approach.lower() != "rmf":
            raise ValueError("Condensing with scheme={0} is only possible in combination with "
                             "approach='RMF'! Given approach={1}".format(scheme, approach))
        if (not hasattr(type(molecule), "__iter__") or len(molecule) != 3 or not
                np.all([isinstance(mol, Molecule) for mol in molecule])):
            raise ValueError("Condensing with scheme={0} needs 3 molecules!".format(scheme))
        # get populations
        pops = [getattr(mol, scheme + "_charges") for mol in molecule]
        if np.any([isinstance(pop, type(None)) for pop in pops]):
            raise ValueError("Condensing scheme={0} is not possible, because attribute {1}_charges "
                             "of molecule instances is 'None'.".format(scheme, scheme.lower()))
        # make dictionary of populations
        dict_pops = dict([(sum(m.mo.nelectrons), m.numbers - pop) for m, pop in zip(molecule, pops)])
        return dict_pops

    # case of condensing the density using denspart
    try:
        # check whether molecules have the same coordinates
        get_matching_attr(molecule, "coordinates", 1.e-4)
        same_coordinates = True
    except ValueError:
        if approach.lower() == "fmr":
            raise ValueError("When geometries of molecules are different, only approach='RMF' "
                             "is possible! Given approach={0}".format(approach.upper()))
        same_coordinates = False

    # find reference molecule
    if isinstance(molecule, Molecule):
        mol0 = molecule
    elif np.all([isinstance(mol, Molecule) for mol in molecule]):
        if len(molecule) != 3:
            raise ValueError("Condensing within FD approach, currently works for "
                             "only 3 molecules! Given {0} molecules.".format(len(molecule)))
        # reference molecule is the middle molecule (for 3 molecules)
        dict_mols = {sum(mol.mo.nelectrons): mol for mol in molecule}
        mol0 = dict_mols.pop(sorted(dict_mols.keys())[1])
    else:
        raise ValueError("Argument molecule not recognized!")

    # check and generate partitioning class & proatomdb
    wpart = wpart_schemes[scheme]
    # make proatom database
    if scheme.lower() not in ["mbis", "b"]:
        if "proatomdb" not in kwargs.keys() or kwargs["proatomdb"] is None:
            proatomdb = ProAtomDB.from_refatoms(mol0.numbers)
            kwargs["proatomdb"] = proatomdb

    # check or generate molecular grid
    grid = get_molecular_grid(molecule, kwargs.pop("grid", None))
    # compute dictionary of number of electron and density
    dict_dens = get_dict_density(molecule, grid.points)

    # compute population of reference molecule
    part0 = wpart(mol0.coordinates, mol0.numbers, mol0.pseudo_numbers, grid,
                  dict_dens[sum(mol0.mo.nelectrons)], **kwargs)
    part0.do_all()
    # record population of reference system
    dict_pops = dict([(sum(mol0.mo.nelectrons), part0["populations"])])
    del dict_dens[sum(mol0.mo.nelectrons)]

    # compute and record populations given grid in a dictionary
    for nelec, dens in dict_dens.iteritems():

        if approach.lower() == "fmr":
            # fragment of molecular response
            pops = condense_to_atoms(dens, part0)
        elif approach.lower() == "rmf":
            # response of molecular fragment
            if not same_coordinates:
                mol = dict_mols[nelec]
                grid = get_molecular_grid(molecule, None)
                dens = molecule.compute_density(grid.points, "ab", None)
            else:
                mol = mol0
            parts = wpart(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                          grid, dens, **kwargs)

            parts.do_all()
            pops = parts["populations"]
        else:
            raise ValueError("Condensing approach {0} is not recognized!".format(approach))
        # Store number of electron and populations in a dictionary
        dict_pops[nelec] = pops
    return dict_pops
