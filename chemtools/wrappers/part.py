import numpy as np
from denspart.mbis import partition
from chemtools.wrappers.molecule import Molecule
from chemtools.wrappers.grid import MolecularGrid
__all__ = ['DensPart']


class DensPart(object):
    """Density-based Atoms-in-Molecules Partitioning Class."""

    available_schemes = ["mbis"]

    def __init__(self, coordinates, numbers, pseudo_numbers,
                 density, grid, scheme):
        if scheme.lower() not in DensPart.available_schemes:
            raise NotImplementedError("MBIS is currently the only supported scheme.")
        self.coordinates = coordinates
        self.numbers = numbers
        self.pseudo_numbers = pseudo_numbers
        self.density = density
        self.grid = grid
        self.part = partition(self.numbers,
                              self.coordinates,
                              self.grid,
                              self.density)
        self.charges = self.part[0].charges

    @classmethod
    def from_molecule(cls, mol, scheme, grid=None, **kwargs):
        if grid is None:
            grid = MolecularGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                                 specs="fine", rotate=False, k=3)
        else:
            check_molecule_grid(mol, grid)
        # compute molecular electron density
        dens = mol.compute_density(grid.points)
        if mol.pseudo_numbers is None:
            pseudo_numbers = mol.numbers
        else:
            pseudo_numbers = mol.pseudo_numbers
        return cls(mol.coordinates, mol.numbers, pseudo_numbers, dens,
                   grid, scheme, **kwargs)

    @classmethod
    def from_file(cls, fname, scheme, grid=None, **kwargs):
        mol = Molecule.from_file(fname)
        return cls.from_molecule(mol, scheme, grid=grid, **kwargs)

    def condense_to_atoms(self, value, w_power=1):
        # promodel, atomic_grids = self.part
        # condensed = np.zeros(promodel.natom)
        # for index in range(promodel.natom):
        #     at_grid = atomic_grids[index]
        #     at_weight = atomic_grids[index].weights
        #     local_prop = promodel.compute_proatom(index, value)
        #     condensed[index] = at_grid.integrate(at_weight**w_power * local_prop)
        # return condensed
        raise NotImplementedError

    def condense_to_fragments(self, value, fragments=None, w_power=1):
        promodel, atomic_grids = self.part
        if fragments is None:
            fragments = [[index] for index in range(promodel.natom)]
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
                weight += atomic_grids[item]
                # weight += self.part.cache.load("at_weights", item)
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
    if not np.max(abs(grid.centers - mol.coordinates)) < 1.e-6:
        raise ValueError("Argument molecule & grid should have the same coordinates/centers.")
    if not np.max(abs(grid.numbers - mol.numbers)) < 1.e-6:
        raise ValueError("Arguments molecule & grid should have the same numbers.")
    if not np.max(abs(grid.pseudo_numbers - mol.pseudo_numbers)) < 1.e-6:
        raise ValueError("Arguments molecule & grid should have the same pseudo_numbers.")