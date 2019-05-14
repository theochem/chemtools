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
# pragma pylint: disable=wrong-import-position
"""Module for (non)bonding interaction analysis of Quantum Chemistry Output Files."""


import numpy as np

from chemtools.wrappers.molecule import Molecule
from chemtools.denstools.densbased import DensGradTool
from chemtools.utils.utils import doc_inherit
from chemtools.utils.cube import CubeGen
from chemtools.outputs.plot import plot_scatter
from chemtools.outputs.vmd import print_vmd_script_nci, print_vmd_script_isosurface

from numpy.ma import masked_less


__all__ = ['NCI', 'ELF', 'LOL']


class BaseInteraction(object):
    """Base class for (non)bonding interactions indicators."""

    @classmethod
    def from_file(cls, fname, spin='ab', index=None, grid=None):
        """Initialize class using wave-function file.

        Parameters
        ----------
        fname : str
            A string representing the path to a molecule's filename.
        spin : str, optional
            The type of occupied spin orbitals; options are 'a', 'b' & 'ab'.
        index : int or Sequence of int, optional
            Sequence of integers representing the index of spin orbitals.
            If None, all occupied spin orbitals are included.
        grid : instance of `Grid`, optional
            Grid used for calculating and visualizing the property values.
            If None, a cubic grid is constructed from molecule with spacing=0.1 & threshold=2.0.
        """
        molecule = Molecule.from_file(fname)
        return cls.from_molecule(molecule, spin=spin, index=index, grid=grid)

    @classmethod
    def from_molecule(cls, molecule, spin='ab', index=None, grid=None):
        """Initialize class from ``Molecule`` object.

        Parameters
        ----------
        molecule : instance of `Molecule` class.
            Instance of `Molecular` class.
        spin : str, optional
            The type of occupied spin orbitals; options are 'a', 'b' & 'ab'.
        index : int or Sequence of int, optional
            Sequence of integers representing the index of spin orbitals.
            If None, all occupied spin orbitals are included.
        grid : instance of `Grid`, optional
            Grid used for calculating and visualizing the property values.
            If None, a cubic grid is constructed from molecule with spacing=0.1 & threshold=2.0.
        """
        pass

    @staticmethod
    def _check_grid(molecule, grid):
        if grid is None:
            grid = CubeGen.from_molecule(molecule.numbers, molecule.pseudo_numbers,
                                         molecule.coordinates, spacing=0.1, threshold=5.0)
        elif not hasattr(grid, 'points'):
            raise ValueError('Argument grid should have "points" attribute!')

        return grid

    @staticmethod
    def _transform(ratio, transformation, k):
        if transformation == 'original':
            return 1.0 / (1.0 + ratio**k)
        elif transformation == 'hyperbolic':
            return NotImplementedError
        elif transformation == 'inverse':
            return NotImplementedError
        else:
            raise ValueError('Argument trans={0} not recognized!'.format(transformation))


class NCI(BaseInteraction):
    """Non-Covalent Interactions (NCI) Class."""

    def __init__(self, density, rdgradient, grid, hessian=None):
        """Initialize class using density, reduced density gradient and `CubeGen` instance.

        Parameters
        ----------
        density : np.array
            Density evaluated on grid points of `cube`.
        rdgradient : np.array
            Reduced density gradient evaluated on grid points of `cube`
        grid : instance of `CubeGen`, optional
            Cubic grid used for calculating and visualizing the NCI.
            If None, it is constructed from molecule with spacing=0.1 and threshold=2.0
        hessian : np.array, optional
            Hessian of density evaluated on grid points of `cube`. This is a array with shape
            (n, 6) where n is the number of grid points of `cube`.
        """
        if density.shape != (len(grid.points),):
            raise ValueError('Shape of density argument {0} does not match '
                             'expected ({1},) shape.'.format(density.shape, len(grid.points)))
        if rdgradient.shape != (len(grid.points),):
            raise ValueError('Shape of rdgradient argument {0} does not '
                             'match expected ({1},) shape.'.format(density.shape, len(grid.points)))

        if hessian is not None:
            if hessian.shape != (len(grid.points), 6):
                raise ValueError('Shape of hessian argument {0} does not match expected ({1}, 6)'
                                 ' shape.'.format(hessian.shape, len(grid.points)))

            # convert the (n, 6) shape to (n, 3, 3) to calculate eigenvalues.
            hestri = np.zeros((len(grid.points), 3, 3))
            tmp = np.zeros((3, 3))
            for i in range(0, len(grid.points)):
                tmp[np.triu_indices(3)] = hessian[i, :]
                hestri[i, :] = tmp

            # compute hessian and its eigenvalues on cubic grid
            eigvalues = np.linalg.eigvalsh(hestri, UPLO='U')

            # use sign of second eigenvalue to distinguish interaction types
            sdens = np.sign(eigvalues[:, 1]) * density

            self._signed_density = sdens
            self._eigvalues = eigvalues

        else:
            self._signed_density = None
            self._eigvalues = None

        self._density = density
        self._rdgrad = rdgradient
        self._grid = grid

    @classmethod
    @doc_inherit(BaseInteraction, 'from_molecule')
    def from_molecule(cls, molecule, spin='ab', index=None, grid=None):
        # generate or check cubic grid
        grid = BaseInteraction._check_grid(molecule, grid)
        # compute density, gradient & hessian on cubic grid
        dens = molecule.compute_density(grid.points, spin=spin, index=index)
        grad = molecule.compute_gradient(grid.points, spin=spin, index=index)
        hess = molecule.compute_hessian(grid.points, spin=spin, index=index)
        # compute reduced gradient
        rdgrad = DensGradTool(dens, grad).reduced_density_gradient
        return cls(dens, rdgrad, grid, hessian=hess)

    @property
    def signed_density(self):
        r"""Signed electron density.

        Electron density :math:`\rho\left(\mathbf{r}\right)` evaluated on a grid,
        signed by the second eigenvalue of the Hessian at that point, i.e.
        :math:`\text{sgn}\left(\lambda_2\right) \times \rho\left(\mathbf{r}\right)`.
        """
        return self._signed_density

    @property
    def eigvalues(self):
        r"""Eigenvalues of Hessian."""
        return self._eigvalues

    def plot(self, fname, color='b'):
        r"""Plot reduced density gradient.

        Reduced density gradient vs.
        :math:`\text{sgn}\left(\lambda_2\right) \times \rho\left(\mathbf{r}\right)`.

        Parameters
        ----------
        fname : str
            A string representing the path to a fname for storing the plot.
            If the given fname does not have a proper extension, the 'png' format is used
            by default, i.e. plot is saved as fname.png.
        color : str, optional
            Color of plot. To customize color, see http://matplotlib.org/users/colors.html

        """
        # scatter plot
        kwargs = {'color': color,
                  'xlim': (-0.2, 0.2),
                  'ylim': (0., 2.),
                  'xlabel': r'sgn$\mathbf{(\lambda_2)}$ $\times$ $\mathbf{\rho(r)}$ (a.u)',
                  'ylabel': 'Reduced Density Gradient'}
        plot_scatter(self._signed_density, self._rdgrad, fname, **kwargs)

    def generate_scripts(self, fname, isosurf=0.50, denscut=0.05):
        r"""Generate cube files and VMD script to visualize non-covalent interactions (NCI).

        Generate density and reduced density gradient cube files, as well as a VMD (Visual
        Molecular Dynamics) script to visualize non-covalent interactions (NCI).

        Parameters
        ----------
        fname : str
            Name of generated cube files and vmd script.
        isosurf : float, optional
            Value of reduced density gradient (RDG) iso-surface used in VMD script.
        denscut : float, optional
            Density cutoff used in creating reduced density gradient cube file.
            Similar to NCIPlot program, reduced density gradient of points with
            density > denscut will be set to 100.0 to display reduced density gradient
            iso-surface subject to the constraint of low density.
            To visualize all reduced density gradient iso-surfaces, disregarding of the
            corresponding density value, set this argument equal to infity using `float('inf')`.
        Note
        ----
        The generated cube files and script imitate the NCIPlot software version 1.0.
        """
        if not isinstance(self._grid, CubeGen):
            raise ValueError("Only possible if argument grid is a cubic grid.")
        # similar to NCIPlot program, reduced density gradient of points with
        # density > cutoff will be set to 100.0 before generating cube file to
        # display reduced density gradient iso-surface subject to the constraint
        # of low density, i.e. density < denscut.
        cutrdg = np.array(self._rdgrad, copy=True)
        cutrdg[abs(self._density) > denscut] = 100.0

        # similar to NCIPlot program, sign(hessian second eigenvalue)*density is
        # multiplied by 100.0 before generating cube file used for coloring the
        # reduced density gradient iso-surface.
        if self._signed_density is not None:
            dens = 100.0 * self._signed_density
        else:
            dens = 100.0 * self._density

        # name of output files
        densfile = fname + '-dens.cube'    # density cube file
        rdgfile = fname + '-grad.cube'     # reduced density gradient cube file
        vmdfile = fname + '.vmd'           # vmd script file
        # dump density & reduced density gradient cube files
        self._grid.dump_cube(densfile, dens)
        self._grid.dump_cube(rdgfile, cutrdg)
        # write VMD scripts
        print_vmd_script_nci(vmdfile, densfile, rdgfile, isosurf, denscut * 100.0)


class ELF(BaseInteraction):
    r"""Electron Localization Function (ELF) introduced by Becke and Edgecombe.

    .. math::
       \text{ELF} (\mathbf{r}) =
            \frac{1}{\left( 1 + \left(\frac{D_{\sigma}(\mathbf{r})}
            {D_{\sigma}^0 (\mathbf{r})} \right)^2\right)}

    with XXX, XXX, and positive definite kinetic energy density defined as, respectively,

    .. math::
        D_{\sigma} (\mathbf{r}) &= \tau_{\sigma} (\mathbf{r}) -
           \frac{1}{4} \frac{(\nabla \rho_{\sigma})^2}{\rho_{\sigma}}

       D_{\sigma}^0 (\mathbf{r}) &=
          \frac{3}{5} (6 \pi^2)^{2/3} \rho_{\sigma}^{5/3} (\mathbf{r})

       \tau_{\sigma} (\mathbf{r}) &=
             \sum_i^{\sigma} \lvert \nabla \phi_i (\mathbf{r}) \rvert^2
    """

    def __init__(self, dens, grad, kin, grid=None, trans='original', k=2):
        r"""Initialize class.

        Parameters
        ----------
        dens : np.ndarray
            Electron density evaluated on a set of grid points.
        grad : np.ndarray
            Gradient vector of electron density evaluated on a set of grid points.
        ked : np.ndarray
            Kinetic energy density evaluated on a set of grid points.
        grid : instance of `Grid`, optional
            Grid used for calculating and visualizing the property values.
            If None, a cubic grid is constructed from molecule with spacing=0.1 & threshold=2.0.
        trans : str, optional
            The approach for transforming the ELF ration to obtain ELF value.
        k : float, optional
            Value of parameter :math:`k`.

        """
        # if dens.shape != (grid.npoints,):
        #     raise ValueError("Arguments dens should have the same size as grid.npoints!")
        # if grad.shape != (grid.npoints, 3):
        #     raise ValueError("Arguments grad should have the same size as grid.npoints!")
        # if kin.shape != (grid.shape,):
        #     raise ValueError("Arguments kin should have the same size as grid.npoints!")
        self._grid = grid
        self._denstool = DensGradTool(dens, grad)
        # compute elf ratio & apply transformation
        self._ratio = kin - self._denstool.ked_weizsacker
        self._ratio /= masked_less(self._denstool.ked_thomas_fermi, 1.0e-30)
        self._value = self._transform(self._ratio, trans, k)
        # assign basins and topology attributes
        self._basins = None
        self._topology = None

    @classmethod
    @doc_inherit(BaseInteraction, 'from_molecule')
    def from_molecule(cls, molecule, spin='ab', index=None, grid=None, trans='original', k=2):
        # generate cubic grid or check grid
        grid = BaseInteraction._check_grid(molecule, grid)
        # compute density, gradient & kinetic energy density on grid
        dens = molecule.compute_density(grid.points, spin=spin, index=index)
        grad = molecule.compute_gradient(grid.points, spin=spin, index=index)
        kin = molecule.compute_kinetic_energy_density(grid.points, spin=spin, index=index)
        return cls(dens, grad, kin, grid, trans, k)

    @property
    def density(self):
        r"""Electron density :math:`\rho\left(\mathbf{r}\right)` evaluated on grid points."""
        return self._denstool.density

    @property
    def ratio(self):
        r"""The ELF ratio evaluated on the grid points."""
        return self._ratio

    @property
    def value(self):
        r"""The :math:`ELF(\mathbf{r})` evaluated on grid points."""
        return self._value

    @property
    def basins(self):
        """Electron localization function basin values of grid points."""
        if not self._basins:
            self._basins = self._compute_basins()
        return self._basins

    @property
    def topology(self):
        """Electron localization function topology."""
        if not self._topology:
            self._topology = self._compute_topology()
        return self._topology

    def _compute_basins(self):
        raise NotImplementedError

    def _compute_topology(self):
        raise NotImplementedError

    def generate_scripts(self, fname, isosurf=0.8):
        """Generate VMD scripts & cube file to visualize ELF iso-surface.

        Parameters
        ----------
        fname : str
            A string representing the path to a fname of generated files.
            The VMD script and cube file will be name fname.vmd and fname-elf.cube, respectively.
        isosurf : float
            Value of ELF iso-surface used in VMD script.

        """
        if not isinstance(self._grid, CubeGen):
            raise ValueError("Only possible if argument grid is a cubic grid.")
        # dump ELF cube file & generate vmd script
        self._grid.dump_cube(fname + '-elf.cube', self.value)
        print_vmd_script_isosurface(fname + '.vmd', fname + '-elf.cube', isosurf=isosurf)


class LOL(BaseInteraction):
    r"""Localized orbital Locator (LOL) Class."""

    def __init__(self, dens, grad, kin, grid, trans='original', k=2):
        self._tool = DensGradTool(dens, grad)
        self._dens = dens
        self._grid = grid
        self._ratio = self._tool.ked_thomas_fermi / kin
        self._value = self._transform(self._ratio, trans, k)

    @classmethod
    @doc_inherit(BaseInteraction, 'from_molecule')
    def from_molecule(cls, molecule, spin='ab', index=None, grid=None, trans='original', k=2):
        # generate cubic grid or check grid
        grid = BaseInteraction._check_grid(molecule, grid)
        # compute density, gradient & kinetic energy density on grid
        dens = molecule.compute_density(grid.points, spin=spin, index=index)
        grad = molecule.compute_gradient(grid.points, spin=spin, index=index)
        kin = molecule.compute_kinetic_energy_density(grid.points, spin=spin, index=index)
        return cls(dens, grad, kin, grid)

    @property
    def ratio(self):
        r"""The LOL ratio evaluated on the grid points."""
        return self._ratio

    @property
    def value(self):
        r"""The :math:`LOL(\mathbf{r})` evaluated on grid points."""
        return self._value

    def generate_scripts(self, fname, isosurf=0.5):
        """Generate VMD scripts & cube file to visualize LOL iso-surface.

        Parameters
        ----------
        fname : str
            A string representing the path to a fname of generated files.
            The VMD script and cube file will be name fname.vmd and fname-lol.cube, respectively.
        isosurf : float
            Value of LOL iso-surface used in VMD script.
        """
        if not isinstance(self._grid, CubeGen):
            raise ValueError("Only possible if argument grid is a cubic grid.")
        # dump LOL cube files
        self._grid.dump_cube(fname + '-lol.cube', self.value)
        # write VMD script for visualization
        print_vmd_script_isosurface(fname + '.vmd', fname + '-lol.cube', isosurf=isosurf)
