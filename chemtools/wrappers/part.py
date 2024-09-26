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
"""Wrapper of Part Module."""
import glob

import numpy as np

# from horton import ProAtomDB
# from horton.scripts.wpart import wpart_schemes

from chemtools.wrappers.molecule import Molecule
from chemtools.wrappers.grid import MolecularGrid

from rhopart import ProAtomDB, ProAtomRecord
from rhopart import VarHirshfeld, DirectedAlphaDivergence, Hirshfeld, HirshfeldI

from grid.onedgrid import GaussChebyshev
from grid.rtransform import BeckeRTransform

from denspart.mbis import MBISProModel
from denspart.vh import optimize_reduce_pro_model
from denspart.properties import compute_radial_moments, compute_multipole_moments, safe_ratio



__all__ = ['DensPart']


class DensPart(object):
    """Density-based Atoms-in-Molecules Partitioning Class."""

    def __init__(self, coordinates, numbers, pseudo_numbers, density, grid, scheme, proatomdb=None, **kwargs):
        """Initialize class.

        Parameters
        ----------
        coordinates : np.ndarray, shape=(M, 3)
            Cartesian coordinates of `M` atoms in the molecule.
        numbers : np.ndarray, shape=(M,)
            Atomic number of `M` atoms in the molecule.
        pseudo_numbers : np.ndarray, shape=(M,)
            Pseudo-number of `M` atoms in the molecule.
        density : np.ndarray, shape=(N,)
            Total density to be partitioned.
        grid : BeckeMolGrid
            Instance of MoleGrid numerical integration grid.
        scheme : str
            Type of atoms-in-molecule partitioning scheme.

        """
        if scheme == "h":
            if not proatomdb:
                proatomdb = ProAtomDB.from_refatoms(numbers)
            part = Hirshfeld(coordinates,
                                  numbers,
                                  pseudo_numbers,
                                  grid,
                                  density,
                                  proatomdb,
                                  lmax=3)
            part.run()
            self.part = part
            self.charges = self.part.charges
            self.at_weights = part.weights
            self.populations = part.populations

        elif scheme == "hi":
            if not proatomdb:
                proatomdb = ProAtomDB.from_refatoms(numbers)
            part = HirshfeldI(coordinates,
                                  numbers,
                                  pseudo_numbers,
                                  grid,
                                  density,
                                  proatomdb,
                                  lmax=3)
            part.run()
            self.part = part
            self.charges = part.charges
            self.at_weights = part.weights
            self.populations = part.populations


        elif scheme == "mbis":
            pro_model_init = MBISProModel.from_geometry(numbers, coordinates)
            pro_model, localgrids = optimize_reduce_pro_model(
                pro_model_init,
                grid,
                density)
            print("Promodel")
            pro_model.pprint()

            print("Computing additional properties")
            results = pro_model.to_dict()
            # Compute atomic weights
            at_weights = np.zeros((numbers.shape[0], grid.points.shape[0]))
            # pro = pro_model.compute_density(grid, localgrids)
            pro = pro_model.compute_density(grid)
            for iatom, atcoord in enumerate(pro_model.atcoords):
                print("Atom index = ", iatom)
                pro_atom = pro_model.compute_proatom(iatom, grid.points)
                # ratio (copied from denspart properties module) is the ratio of molecular and promolecular density
                # ratio = safe_ratio(density[localgrid.indices], pro[localgrid.indices])
                # atomic weight
                atweight = safe_ratio(pro_atom, pro)
                # check: compute atomic population to make sure it matches expected values (from Denspart)
                print(" Computed Charge = ", numbers[iatom] - grid.integrate(
                    atweight * density))
                print(" Expected Charge = ", pro_model.charges[iatom])
                print("")
                print(atweight)
                at_weights[iatom] = atweight

            self.part = None
            self.charges = pro_model.charges
            self.at_weights = at_weights
            self.populations = numbers - pro_model.charges




        self.grid = grid
        self.density = density
        self.coordinates = coordinates
        self.numbers = numbers
        self.pseudo_numbers = pseudo_numbers



    @classmethod
    def from_molecule(cls, mol, scheme=None, grid=None, spin="ab", **kwargs):
        """Initialize class given a Molecule instance.

        Parameters
        ----------
        mol : Molecule
            Instance of Molecule class.
        scheme : str
            Type of atoms-in-molecule partitioning scheme.
        grid : MolecularGrid
            Instance of MolecularGrid numerical integration grid.
        spin : str
           Type of occupied spin orbitals; choose either "a" (for alpha), "b" (for beta),
           and "ab" (for alpha + beta).

        """
        if grid is None:
            oned = GaussChebyshev(60)
            rgrid = BeckeRTransform(1e-5, 1).transform_1d_grid(oned)
            grid = MolecularGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                                 specs=[rgrid, 590], rotate=False)
        else:
            check_molecule_grid(mol, grid)
        # compute molecular electron density
        dens = mol.compute_density(grid.points, spin=spin)
        if mol.pesudo_numbers is None:
            pesudo_numbers = mol.numbers
        else:
            pesudo_numbers = mol.pesudo_numbers
        return cls(mol.coordinates, mol.numbers, pesudo_numbers, dens, grid, scheme, **kwargs)

    @classmethod
    def from_file(cls, fname, scheme=None, grid=None, spin="ab", **kwargs):
        """Initialize class given a wave-function file.

        Parameters
        ----------
        fname : str
            Path to molecule's files.
        scheme : str
            Type of atoms-in-molecule partitioning scheme.
        grid : MolecularGrid
            Instance of MolecularGrid integration grid.
        spin : str
           Type of occupied spin orbitals; choose either "a" (for alpha), "b" (for beta),
           and "ab" (for alpha + beta).

        """
        mol = Molecule.from_file(fname)
        return cls.from_molecule(mol, scheme=scheme, grid=grid, spin=spin, **kwargs)

    def condense_to_atoms(self, value, w_power=1, local=False, coords=None, r=5):
        if not local:
            at_weights = self.at_weights
            prop = (at_weights**w_power) * value[None, :]
            condensed = np.array([self.grid.integrate(prop[i]) for i in range(len(self.numbers))])
        else:
            at_weights_all = self.at_weights
            condensed = []
            for a in range(self.numbser):
                atgrid = self.grid._grid.get_localgrid(coords[a], r)
                at_weights = at_weights_all[a, atgrid.indices]
                prop = (at_weights**w_power) * value[atgrid.indices]
                condensed.append(atgrid.integrate(prop))
            condensed = np.array(condensed)
        return condensed

    def condense_to_fragments(self, value, fragments=None, w_power=1):
        if fragments is None:
            fragments = [[index] for index in range(len(self.numbers))]
        else:
            # check fragments
            segments = sorted([item for frag in fragments for item in frag])
            segments = np.array(segments)
            if segments.size != self.numbers.size:
                raise ValueError("Items in Fragments should uniquely represent all atoms.")
        condensed = np.zeros(len(fragments))
        for index, frag in enumerate(fragments):
            weight = np.zeros(self.grid.points.shape[0])
            for item in frag:
                weight += self.at_weights[item]
            share = self.grid.integrate(weight ** w_power, value)
            condensed[index] = share
        return condensed


def check_molecule_grid(mol, grid):
    """Check that molecule and grid represent the same system.

    Parameters
    ----------
    mol : Molecule
        Instance of Molecule class.
    grid : MolecularGrid
        Instance of MolecularGrid numerical integration grid.

    """
    if not np.max(abs(grid.atcoords - mol.coordinates)) < 1.e-6:
        raise ValueError("Argument molecule & grid should have the same coordinates/centers.")
    # if not np.max(abs(grid.numbers - mol.numbers)) < 1.e-6:
    #     raise ValueError("Arguments molecule & grid should have the same numbers.")
    # if not np.max(abs(grid.pseudo_numbers - mol.pseudo_numbers)) < 1.e-6:
    #     raise ValueError("Arguments molecule & grid should have the same pseudo_numbers.")
