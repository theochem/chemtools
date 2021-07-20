import numpy as np
from denspart.mbis import partition
from chemtools.wrappers.molecule import Molecule
from chemtools.wrappers.grid import MolecularGrid
__all__ = ['DensPart']


class DensPart(object):
    """Density-based Atoms-in-Molecules Partitioning Class."""

    available_schemes = ["mbis"]

    def __init__(self, coordinates, numbers, pseudo_numbers,
                 density, grid, scheme="mbis", **kwargs):
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
                              self.density
                              **kwargs)

    @classmethod
    def from_molecule(cls, mol, scheme=None, grid=None, **kwargs):
        if grid is None:
            grid = MolecularGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                                 specs="fine", rotate=False, k=3)
        else:
            check_molecule_grid(mol, grid)
        # compute molecular electron density
        dens = mol.compute_density(grid.points)
        if mol.pseudo_numbers is None:
            mol.pseudo_numbers = mol.numbers
        else:
            mol.pseudo_numbers = mol.pseudo_numbers
        return cls(mol.coordinates, mol.numbers, mol.pseudo_numbers, dens, grid, scheme, **kwargs)

    @classmethod
    def from_file(cls, fname, scheme=None, grid=None, **kwargs):
        mol = Molecule.from_file(fname)
        return cls.from_molecule(mol, scheme=scheme, grid=grid, **kwargs)

    def condense_to_atoms(self, value, w_power=1):
        raise NotImplementedError

    def condense_to_fragments(self, value, fragments=None, w_power=1):
        raise NotImplementedError

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