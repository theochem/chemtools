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


import sys
import numpy as np

from chemtools.denstools.densbased import DensGradTool
from chemtools.utils.utils import doc_inherit
from chemtools.utils.cube import UniformGrid
from chemtools.outputs.plot import plot_scatter
from chemtools.outputs.vmd import print_vmd_script_nci, print_vmd_script_isosurface

from numpy.ma import masked_less

if sys.version_info.major == 2:
    from chemtools.wrappers2.molecule import Molecule


class BaseInteraction(object):
    """Base class for (non)bonding interactions indicators."""

    @classmethod
    def from_file(cls, fname, spin='ab', index=None, grid=None):
        """Initialize class using wave-function file.

        Parameters
        ----------
        fname : str
            A string representing the path to a molecule's fname.
        spin : str, optional
            The type of occupied spin orbitals; options are 'a', 'b' & 'ab'.
        index : int or Sequence of int, optional
            Sequence of integers representing the index of spin orbitals.
            If None, all occupied spin orbitals are included.
        grid : instance of `Grid`, optional
            Grid used for calculating and visualizing the property values.
            If None, a cubic grid is constructed from molecule with spacing=0.1 & extension=2.0.
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
            If None, a cubic grid is constructed from molecule with spacing=0.1 & extension=2.0.
        """
        pass

    @staticmethod
    def _check_grid(molecule, grid):
        if grid is None:
            grid = UniformGrid.from_molecule(molecule, spacing=0.1, extension=2.0)
        elif not hasattr(grid, 'points'):
            raise ValueError('Argument grid should have "points" attribute!')

        return grid

    @staticmethod
    def _transform(ratio, trans, trans_k, trans_a):
        if trans == 'rational':
            return 1.0 / (1.0 + trans_a * ratio ** trans_k)
        elif trans == 'hyperbolic':
            return 0.5 * (1 + np.tanh(trans_a * (ratio ** -trans_k - ratio ** trans_k)))
        elif trans == 'inverse_rational':
            return 1.0 - 1.0 / (1.0 + trans_a * ratio ** trans_k)
        elif trans == 'inverse_hyperbolic':
            return 0.5 * (1 + np.tanh(-trans_a * (ratio ** -trans_k - ratio ** trans_k)))
        else:
            raise ValueError('Argument trans={0} not recognized!'.format(trans))


class NCI(BaseInteraction):
    """Non-Covalent Interactions (NCI) Class."""

    def __init__(self, density, rdgradient, grid, hessian=None):
        """Initialize class using density, reduced density gradient and `UniformGrid` instance.

        Parameters
        ----------
        density : np.array
            Density evaluated on grid points of `cube`.
        rdgradient : np.array
            Reduced density gradient evaluated on grid points of `cube`
        grid : instance of `UniformGrid`, optional
            Cubic grid used for calculating and visualizing the NCI.
            If None, it is constructed from molecule with spacing=0.1 and extension=2.0.
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
            if hessian.shape != (len(grid.points), 3, 3):
                raise ValueError("Shape of hessian argument {0} does not match expected "
                                 "({1}, 3, 3) shape!".format(hessian.shape, len(grid.points)))

            # compute hessian eigenvalues on cubic grid
            eigvalues = np.linalg.eigvalsh(hessian, UPLO='U')

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

    def generate_plot(self, fname, color='b', denslim=(-0.2, 0.2), rdglim=(0., 2.)):
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
        denslim: tuple, optional
            The minimum and maximum of the (signed) density in the plot.
        rdglim: tuple, optional
            The minimum and maximum of the reduced density gradient in the plot.

        """
        # scatter plot
        kwargs = {'color': color,
                  'xlim': denslim,
                  'ylim': rdglim,
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
            corresponding density value, set this argument equal to infinity using `float('inf')`.
        Note
        ----
        The generated cube files and script imitate the NCIPlot software version 1.0.
        """
        if not isinstance(self._grid, UniformGrid):
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
        self._grid.generate_cube(densfile, dens)
        self._grid.generate_cube(rdgfile, cutrdg)
        # write VMD scripts
        print_vmd_script_nci(vmdfile, densfile, rdgfile, isosurf, denscut * 100.0)


class ELF(BaseInteraction):
    r"""Electron Localization Function (ELF) introduced by Becke and Edgecombe.

    The ELF ratio is defined as:

    .. math::
       \zeta_\text{ELF}(\mathbf{r}) = \frac{\tau_\text{PD}(\mathbf{r}) -
       \tau_\text{W}(\mathbf{r})}{\tau_\text{TF}(\mathbf{r})}

    where :math:`\tau_\text{PD}(\mathbf{r})`, :math:`\tau_\text{W}(\mathbf{r})` and
    :math:`\tau_\text{TF}(\mathbf{r})` are positive-definite (Lagrangian), Weizsacker and
    Thomas-Fermi kinetic energy densities defined in :class:`chemtools.toolbox.kinetic.KED`.

    The ELF is computed by transforming the ratio:

    .. math:: \text{ELF}(\mathbf{r}) = f\left(\zeta_\text{ELF}(\mathbf{r})\right)

    where the transformation :math:`f` can be:

    .. math::
       \text{rational  : } \, f(\zeta, k, a) &= \frac{1}{1 + a \, \zeta^k} \\
       \text{hyperbolic: } \, f(\zeta, k, a) &= \tfrac{1}{2}
              \left(1 + \tanh\left(a \left(\zeta^{-k} - \zeta^{k}\right)\right)\right)

    Traditionally, the **'rational'** transformation with :math:`k=2` and :math:`a=1` is used.

    """

    def __init__(self, dens, grad, ked, grid=None, trans='rational', trans_k=2, trans_a=1,
                 denscut=0.0005):
        r"""Initialize class from arrays.

        Parameters
        ----------
        dens : np.ndarray
            Electron density of grid points, :math:`\rho(\mathbf{r})`.
        grad : np.ndarray
            Gradient of electron density of grid points, :math:`\nabla \rho(\mathbf{r})`.
        ked : np.ndarray
            Positive-definite or Lagrangian kinetic energy density of grid
            points; :math:`\tau_\text{PD} (\mathbf{r})` or :math:`G(\mathbf{r})`.
        grid : instance of `Grid`, optional
            Grid used for computation of ELF. Only if this a CubeGrid one can generate the scripts.
        trans : str, optional
            Type of transformation applied to ELF ratio; options are 'rational' or 'hyperbolic'.
        trans_k : float, optional
            Parameter :math:`k` of transformation.
        trans_a : float, optional
            Parameter :math:`a` of transformation.
        denscut : float, optional
            Value of density cut. ELF value of points with density < denscut is set to zero.

        """
        if dens.shape != ked.shape:
            raise ValueError('Arguments dens and ked should have the same shape!')
        if grad.ndim != 2:
            raise ValueError('Argument grad should be a 2d-array!')
        if grad.shape[0] != dens.shape[0]:
            raise ValueError('Argument dens & grad should have the same length!')
        if trans.lower() not in ['rational', 'hyperbolic']:
            raise ValueError('Argument trans should be either "rational" or "hyperbolic".')
        if not trans_k > 0:
            raise ValueError('Argument trans_k should be positive! trans_k={0}'.format(trans_k))
        if not trans_a > 0:
            raise ValueError('Argument trans_a should be positive! trans_a={0}'.format(trans_a))
        self._grid = grid
        self._denstool = DensGradTool(dens, grad)
        # compute elf ratio
        self._ratio = ked - self._denstool.ked_weizsacker
        self._ratio /= masked_less(self._denstool.ked_thomas_fermi, 1.0e-30)
        # compute elf value & set low density points to zero
        self._value = np.asarray(self._transform(self._ratio, trans.lower(), trans_k, trans_a))
        self._value[self._denstool.density < denscut] = 0.

    @classmethod
    def from_molecule(cls, molecule, spin='ab', index=None, grid=None, trans='rational',
                      trans_k=2, trans_a=1, denscut=0.0005):
        """Initialize class from molecule.

        Parameters
        ----------
        molecule : instance of `Molecule` class.
            Instance of `Molecular` class.
        spin : str, optional
            Type of occupied spin orbitals; options are 'a', 'b' & 'ab'.
        index : int or sequence of int, optional
            Sequence of spin orbital indices to use. If None, all occupied spin orbitals are used.
        grid : instance of `Grid`, optional
            Grid used for computation of ELF. Only if this a CubeGrid one can generate the scripts.
            If None, a cubic grid is constructed from molecule with spacing=0.1 & extension=2.0.
        trans : str, optional
            Type of transformation applied to ELF ratio; options are 'rational' or 'hyperbolic'.
        trans_k : float, optional
            Parameter :math:`k` of transformation.
        trans_a : float, optional
            Parameter :math:`a` of transformation.
        denscut : float, optional
            Value of density cut. ELF value of points with density < denscut is set to zero.

        """
        # generate cubic grid or check grid
        grid = BaseInteraction._check_grid(molecule, grid)
        # compute density, gradient & kinetic energy density on grid
        dens = molecule.compute_density(grid.points, spin=spin, index=index)
        grad = molecule.compute_gradient(grid.points, spin=spin, index=index)
        kin = molecule.compute_ked(grid.points, spin=spin, index=index)
        return cls(dens, grad, kin, grid, trans, trans_k, trans_a, denscut)

    @classmethod
    def from_file(cls, fname, spin='ab', index=None, grid=None, trans='rational',
                  trans_k=2, trans_a=1, denscut=0.0005):
        """Initialize class from wave-function file.

        Parameters
        ----------
        fname : str
            Path to a molecule's file.
        spin : str, optional
            Type of occupied spin orbitals; options are 'a', 'b' & 'ab'.
        index : int or sequence of int, optional
            Sequence of spin orbital indices to use. If None, all occupied spin orbitals are used.
        grid : instance of `Grid`, optional
            Grid used for computation of ELF. Only if this a CubeGrid one can generate the scripts.
            If None, a cubic grid is constructed from molecule with spacing=0.1 & extension=2.0.
        trans : str, optional
            Type of transformation applied to ELF ratio; options are 'rational' or 'hyperbolic'.
        trans_k : float, optional
            Parameter :math:`k` of transformation.
        trans_a : float, optional
            Parameter :math:`a` of transformation.
        denscut : float, optional
            Value of density cut. ELF value of points with density < denscut is set to zero.

        """
        molecule = Molecule.from_file(fname)
        return cls.from_molecule(molecule, spin, index, grid, trans, trans_k, trans_a, denscut)

    @property
    def ratio(self):
        r"""The ELF ratio evaluated on grid points."""
        return self._ratio

    @property
    def value(self):
        r"""The :math:`\text{ELF}(\mathbf{r})` evaluated on grid points."""
        return self._value

    def generate_scripts(self, fname, isosurf=0.8):
        """Generate VMD scripts & cube file to visualize ELF iso-surface.

        Parameters
        ----------
        fname : str
            File name used for the generated files.
            The VMD script and cube file will be named fname.vmd and fname-elf.cube, respectively.
        isosurf : float, optional
            Value of ELF iso-surface used in VMD script.

        """
        if not isinstance(self._grid, UniformGrid):
            raise ValueError('Only possible if argument grid is a cubic grid.')
        if self._denstool.density.shape[0] != self._grid.points.shape[0]:
            raise ValueError('Number of grid points should match number of dens values!')
        # dump ELF cube file & generate vmd script
        vmdname = fname + '.vmd'
        cubname = fname + '-elf.cube'
        self._grid.generate_cube(cubname, self.value)
        print_vmd_script_isosurface(vmdname, cubname, isosurf=isosurf, representation='Line')


class LOL(BaseInteraction):
    r"""Localized Orbital Locator (LOL) introduced by Becke and Schmider.

    The LOL ratio is defined as:

    .. math::
       \zeta_\text{LOL}(\mathbf{r}) = \frac{\tau_\text{TF}(\mathbf{r})}{\tau_\text{PD}(\mathbf{r})}

    where :math:`\tau_\text{TF}(\mathbf{r})` and :math:`\tau_\text{PD}(\mathbf{r})` are
    Thomas-Fermi and positive-definite (Lagrangian) kinetic energy densities defined in
    :class:`chemtools.toolbox.kinetic.KED`.

    The LOL is computed by transforming the ratio:

    .. math:: \text{LOL}(\mathbf{r}) = f\left(\zeta_\text{LOL}(\mathbf{r})\right)

    where the transformation :math:`f` can be:

    .. math::
       \text{inverse_rational  : } \, f(\zeta, k, a) &= 1 - \frac{1}{1 + a \, \zeta^k} \\
       \text{inverse_hyperbolic: } \, f(\zeta, k, a) &= \tfrac{1}{2}
              \left(1 + \tanh\left(-a \left(\zeta^{-k} - \zeta^{k}\right)\right)\right)

    Traditionally, the **'inverse_rational'** transformation with :math:`k=1` and :math:`a=1`
    is used.

    """

    def __init__(self, dens, grad, ked, grid=None, trans='inverse_rational', trans_k=1, trans_a=1,
                 denscut=0.0005):
        r"""Initialize class from arrays.

        Parameters
        ----------
        dens : np.ndarray
            Electron density of grid points, :math:`\rho(\mathbf{r})`.
        grad : np.ndarray
            Gradient of electron density of grid points, :math:`\nabla \rho(\mathbf{r})`.
        ked : np.ndarray
            Positive-definite or Lagrangian kinetic energy density of grid
            points; :math:`\tau_\text{PD} (\mathbf{r})` or :math:`G(\mathbf{r})`.
        grid : instance of `Grid`, optional
            Grid used for computation of LOL. Only if this a CubeGrid one can generate the scripts.
        trans : str, optional
            Type of transformation applied to LOL ratio; options are 'inverse_rational' or
            'inverse_hyperbolic'.
        trans_k : float, optional
            Parameter :math:`k` of transformation.
        trans_a : float, optional
            Parameter :math:`a` of transformation.
        denscut : float, optional
            Value of density cut. LOL value of points with density < denscut is set to zero.

        """
        if dens.shape != ked.shape:
            raise ValueError('Arguments dens and ked should have the same shape!')
        if grad.ndim != 2:
            raise ValueError('Argument grad should be a 2d-array!')
        if grad.shape[0] != dens.shape[0]:
            raise ValueError('Argument dens & grad should have the same length!')
        if trans.lower() not in ['inverse_rational', 'inverse_hyperbolic']:
            raise ValueError('Argument trans should be either "inverse_rational" or '
                             '"inverse_hyperbolic".')
        if not trans_k > 0:
            raise ValueError('Argument trans_k should be positive! trans_k={0}'.format(trans_k))
        if not trans_a > 0:
            raise ValueError('Argument trans_a should be positive! trans_a={0}'.format(trans_a))
        self._denstool = DensGradTool(dens, grad)
        self._grid = grid
        # compute elf ratio
        self._ratio = self._denstool.ked_thomas_fermi / masked_less(ked, 1.0e-30)
        # compute elf value & set low density points to zero
        self._value = np.asarray(self._transform(self._ratio, trans.lower(), trans_k, trans_a))
        self._value[self._denstool.density < denscut] = 0

    @classmethod
    def from_molecule(cls, molecule, spin='ab', index=None, grid=None, trans='inverse_rational',
                      trans_k=1, trans_a=1, denscut=0.0005):
        """Initialize class from molecule.

        Parameters
        ----------
        molecule : instance of `Molecule` class.
            Instance of `Molecular` class.
        spin : str, optional
            Type of occupied spin orbitals; options are 'a', 'b' & 'ab'.
        index : int or sequence of int, optional
            Sequence of spin orbital indices to use. If None, all occupied spin orbitals are used.
        grid : instance of `Grid`, optional
            Grid used for computation of LOL. Only if this a CubeGrid one can generate the scripts.
            If None, a cubic grid is constructed from molecule with spacing=0.1 & extension=2.0.
        trans : str, optional
            Type of transformation applied to ELF ratio; options are 'inverse_rational' or
            'inverse_hyperbolic'.
        trans_k : float, optional
            Parameter :math:`k` of transformation.
        trans_a : float, optional
            Parameter :math:`a` of transformation.
        denscut : float, optional
            Value of density cut. LOL value of points with density < denscut is set to zero.

        """
        # generate cubic grid or check grid
        grid = BaseInteraction._check_grid(molecule, grid)
        # compute density, gradient & kinetic energy density on grid
        dens = molecule.compute_density(grid.points, spin=spin, index=index)
        grad = molecule.compute_gradient(grid.points, spin=spin, index=index)
        ked = molecule.compute_ked(grid.points, spin=spin, index=index)
        return cls(dens, grad, ked, grid, trans, trans_k, trans_a, denscut)

    @classmethod
    def from_file(cls, fname, spin='ab', index=None, grid=None, trans='inverse_rational',
                  trans_k=1, trans_a=1, denscut=0.0005):
        """Initialize class from wave-function file.

        Parameters
        ----------
        fname : str
            Path to a molecule's file.
        spin : str, optional
            Type of occupied spin orbitals; options are 'a', 'b' & 'ab'.
        index : int or sequence of int, optional
            Sequence of spin orbital indices to use. If None, all occupied spin orbitals are used.
        grid : instance of `Grid`, optional
            Grid used for computation of LOL. Only if this a CubeGrid one can generate the scripts.
            If None, a cubic grid is constructed from molecule with spacing=0.1 & extension=2.0.
        trans : str, optional
            Type of transformation applied to LOL ratio; options are 'inverse_rational' or
            'inverse_hyperbolic'.
        trans_k : float, optional
            Parameter :math:`k` of transformation.
        trans_a : float, optional
            Parameter :math:`a` of transformation.
        denscut : float, optional
            Value of density cut. LOL value of points with density < denscut is set to zero.

        """
        molecule = Molecule.from_file(fname)
        return cls.from_molecule(molecule, spin, index, grid, trans, trans_k, trans_a, denscut)

    @property
    def ratio(self):
        r"""The LOL ratio evaluated on the grid points."""
        return self._ratio

    @property
    def value(self):
        r"""The :math:`\text{LOL}(\mathbf{r})` evaluated on grid points."""
        return self._value

    def generate_scripts(self, fname, isosurf=0.5):
        """Generate VMD scripts & cube file to visualize LOL iso-surface.

        Parameters
        ----------
        fname : str
            A string representing the path to a fname of generated files.
            The VMD script and cube file will be named fname.vmd and fname-lol.cube, respectively.
        isosurf : float
            Value of LOL iso-surface used in VMD script.

        """
        if not isinstance(self._grid, UniformGrid):
            raise ValueError("Only possible if argument grid is a cubic grid.")
        # dump LOL cube files
        fname_vmd = fname + '.vmd'
        fname_lol = fname + '-lol.cube'
        self._grid.generate_cube(fname_lol, self.value)
        # write VMD script for visualization
        print_vmd_script_isosurface(fname_vmd, fname_lol, isosurf=isosurf, representation='Line')
