# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2023 The ChemTools Development Team
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
"""Module for Interacting quantum atoms(IQA)."""

import logging

import copy as copy
import numpy as np
from numpy.testing import assert_almost_equal

from chemtools.toolbox.utils import check_arg_molecule, \
                                    get_horton_analytical_components, \
                                    get_molecular_grid, get_libxc_xc_density, \
                                    compute_molecular_orbitals_from_ao
from gbasis.wrappers import from_iodata
from gbasis.evals.density import evaluate_density
from gbasis.evals.density import evaluate_general_kinetic_energy_density
from gbasis.evals.density import evaluate_posdef_kinetic_energy_density
from gbasis.evals.density import evaluate_density_gradient
from gbasis.evals.eval import evaluate_basis
from gbasis.evals.eval_deriv import evaluate_deriv_basis
from gbasis.integrals.point_charge import point_charge_integral
from chemtools.utils.cube import UniformGrid
from chemtools.wrappers.grid import MolecularGrid
from chemtools.wrappers.molecule import Molecule
from chemtools.wrappers.part import DensPart

from grid.ngrid import Ngrid
from grid.molgrid import MolGrid

from horton.meanfield import RLibXCLDA, RLibXCGGA, RLibXCHybridGGA, RLibXCMGGA, RLibXCHybridMGGA
from horton.cache import Cache

class IQA(object):
    """Interacting Quantum Atoms (IQA) Class."""

    # def __init__(self, molecule, basis, dm, grid, part=None, molecule_chemtools=None):
    def __init__(self, molecule, basis, dm, grid, part=None, grid_2=None, part_2=None):
        """Initialize class.

        Parameters
        ----------
        molecule : `Molecule`
            Instance of `Molecular` class from Chemtools.
        basis: tuple of gbasis.contraction.GeneralizedContractionShell
            Basis set object used within the `gbasis` module.
        `   GeneralizedContractionShell` corresponds to the `Shell` object within `iodata.basis`.
            A list of strings, each on being either ``"cartesian"`` or  ``"spherical"``
        dm: np.ndarray(nAO, nAO)
            One-electron density matrix in the basis of the Atomic Orbitals.
        grid : `grid`
            Instance of `Grid` class from Chemtools.
        part : `DensePart`, optional
            Instance of `DensePart` class from Chemtools.

        """
        # check basis
        if not basis[0].__class__.__name__ == "IODataShell":
            raise TypeError("basis should have been created with Gbasis")

        # check dm
        if not (isinstance(dm, np.ndarray) and dm.ndim == 2 and dm.dtype == float):
            raise TypeError("One-electron density matrix must be a 2D-array float `dtype`")
        if dm.shape[0] != dm.shape[1]:
            raise ValueError("One-electron density matrix must be a square matrix.")
        if not np.allclose(dm, dm.T):
            raise ValueError("One-electron density matrix must be symmetric.")

        # check grid and part
        if part is not None:
            if not isinstance(part, DensPart) and part.__class__.__name__ not in ['VarHirshfeld', 'HirshfeldI','Hirshfeld']:
                raise TypeError('Argument part should be an instance of DensPart class or VarHirshfeld from rhopart.')
        if part is not None and not (part.numbers == molecule.atnums).all():
            raise ValueError("DensPart molecule different from molecule")
        if not isinstance(grid, UniformGrid) and not isinstance(grid, MolecularGrid):
            raise TypeError(
                "Argument part should be an instance of MolecularGrid or UniformGrid class."
            )

        # Check Grid
        if not isinstance(grid, UniformGrid) and not isinstance(grid, MolecularGrid):
            raise TypeError(
                f"Argument grid should be an instance of MolecularGrid or UniformGrid class. Got {grid.__class__.__name__}"
            )
        # if not (molecule.coordinates == grid.centers).all():
        #     raise ValueError("Molecule and Grid initialized from different molecules")
        # Check optional grid_2 and part_2
        if grid_2:
            if not isinstance(grid_2, UniformGrid) and not isinstance(grid_2, MolecularGrid):
                raise TypeError(
                    "Argument part should be an instance of MolecularGrid or UniformGrid class."
                )
            if not (grid_2.numbers == molecule.atnums).all():
                raise ValueError("Grid_2 molecule different from molecule")
        if part_2 is not None:
            if not isinstance(part, DensPart) and part.__class__.__name__ not in ['VarHirshfeld', 'HirshfeldI','Hirshfeld']:
                raise TypeError('Argument part should be an instance of DensPart class or VarHirshfeld from rhopart.')
        if part_2 is not None and not (part.numbers == molecule.atnums).all():
            raise ValueError("DensPart molecule different from molecule")

        self.molecule = molecule
        self.basis = basis
        self.dm = dm
        self.grid = grid
        self.part = part
        self.grid_2 = grid_2
        self.part_2 = part_2
        dens = evaluate_density(dm, basis, grid.points)
        self.dens = dens

    @classmethod
    def from_molecule(cls, molecule, grid, part=None, scheme=None):
        """Initialize Interacting Quantum Atoms (IQA) class from `Molecules` object.

        Parameters
        ----------
        molecule : `Molecule`
             Instance of `Molecular` class from Chemtools.
        grid : `grid`
            Instance of `Grid` class from Chemtools.
        part : `DensePart`, optional
            Instance of `DensePart` class from Chemtools. If not provided DensPart object will be
            instantiated using scheme.
        scheme : str, optional
            Name of the atomic partition scheme.Default value is None and no atomic decomposition
            of the IQA components is performed. If not provided a part argument corresponding to
            DensPart class object must be passed.

        """
        molecule = check_arg_molecule(molecule)
        # Check restricted closed shell single-determinant wave-function (molecule_iodata)
        if not isinstance(molecule, Molecule):
            raise ValueError("`molecule` must be an instance of class `Molecule` in Chemtools.")
        if molecule._iodata.mo.kind != "restricted":
            raise ValueError("Currently, only 'restricted' wave-functions are supported")
        if int(np.sum(molecule._mo._occs_a)) != int(np.sum(molecule._mo._occs_b)):
            raise ValueError("Code is not tested for open-shell wave-functions.")

        # Check Grid
        if not isinstance(grid, UniformGrid) and not isinstance(grid, MolecularGrid):
            raise TypeError(
                f"Argument grid should be an instance of MolecularGrid or UniformGrid class. Got {grid.__class__.__name__}"
            )
        if not (molecule.coordinates == grid.centers).all():
            raise ValueError("Molecule and Grid initialized from different molecules")

        # Check/Initialize DensPart object
        if part is None and scheme is None:
            print("No atomic partition scheme provided. No atomic decomposition will be performed.")
        elif part is not None:
            if not isinstance(part, DensPart):
                raise TypeError("Argument part should be an instance of DensPart class.")
        elif scheme == 'H':
            part = DensPart.from_molecule(molecule, grid=grid, scheme="h",
                                              local=False)
        elif scheme == 'HI':
            part = DensPart.from_molecule(molecule, grid=grid, scheme="hi",
                                              local=False)
        elif scheme is not None and scheme not in ['H', 'HI']:
            raise NotImplementedError(f"Atomic partition {scheme} not yet available")

        # Initialize gbasis
        basis = from_iodata(molecule._iodata)
        one_rdm = molecule._iodata.one_rdms.get("post_scf", molecule._iodata.one_rdms.get("scf"))
        # Check if dm present
        if one_rdm is None:
            # Check if mo present to construct dm
            if molecule.mo is None:
                raise NotImplementedError("Missing molecule.mo object.Density matrix...")
            else:
                print("Couldn't read Density matrix. Calculating from its components.")
                one_rdm = molecule._mo.compute_dm()

        return cls(molecule, basis, one_rdm, grid, part)

    @classmethod
    def from_file(cls, fname, grid_type=None, scheme=None):
        """Initialize Interacting Quantum Atoms (IQA) class from wave-function file.

        Parameters
        ----------
        fname : str
            Path to molecule's wave-function file.
        grid_type : str, optional
            Type of molecular grid. Options available: Becke, Cube
        scheme : str, optional
            Name of the atomic partition scheme.Default value is None and no atomic decomposition
            of the IQA components is performed.
        """

        molecule = Molecule.from_file(fname)

        # Initialize grid
        # TODO: Add grid type argument to get_molecular_grid
        grid = get_molecular_grid(molecule)

        return cls.from_molecule(molecule, grid, scheme=scheme)

    def iqa(self, dft_exch=None, dft_corr=None):
        """Return the Interacting Quantum Atoms (IQA) components integrated at the evaluated given points

        Parameters
        ----------
        dft_exch: str, optional
            DFT exchange type to obtain Exchange energy density.
        dft_corr: str, optional
            DFT correlation type to obtain Exchange energy density.

        Returns
        -------
        iqa_results: dict
            Dictionary with IQA components.

        """
        molecule = self.molecule
        basis = self.basis
        dm = self.dm
        grid = self.grid
        part = self.part
        rho = self.dens

        # Initialize results dict
        iqa_results = {}

        logging.info('INITIALIZING INTERACTING QUANTUM ATOMS(IQA) CALCULATION')
        iqa_results['nn_total'] = self.nn_iqa()
        iqa_results['en_total'], iqa_results['en_atomic'] = self.en_iqa()
        iqa_results['kin_total'], iqa_results['kin_atomic'], \
            iqa_results['kin_total_posdef'], iqa_results['kin_atomic_posdef']= self.kin_iqa()
        analytical_comp = get_horton_analytical_components(molecule)
        nn_horton = analytical_comp['nn']
        ne_horton = np.trace(dm.dot(analytical_comp['na']))
        kinetic_horton = np.trace(dm.dot(analytical_comp['kin']))
        # todo: test gbasis

        dft_xc_edens = {}
        # assuming dft_corr and dft_exch specified together
        if dft_corr and dft_exch:
            dft_xc_edens.update(get_libxc_xc_density(molecule, grid, analytical_comp, dft_exch, dft_corr))
            iqa_results['x_hf_total'], iqa_results['coul_total'], iqa_results['x_hf_atomic'], \
            iqa_results['coul_atomic'] = self.ee_iqa_hf()
            iqa_results['c_total'], iqa_results['c_atomic'] = self.ee_iqa_dft(rho, dft_xc_edens["edens_c"])
            iqa_results['x_total'], iqa_results['x_atomic'] = self.ee_iqa_dft(rho, dft_xc_edens["edens_x"])
            if 'coeff_mix' in dft_xc_edens.keys() and dft_xc_edens['coeff_mix']:
                if dft_xc_edens['coeff_mix']:
                    # Scaling HF exchange
                    iqa_results['x_hf_total'] = iqa_results['x_hf_total'] * dft_xc_edens[
                        'coeff_mix']
            else:
                iqa_results.pop('x_hf_total')
                iqa_results.pop('x_hf_atomic')

        elif dft_exch:
            dft_xc_edens.update(
                get_libxc_xc_density(molecule, grid, analytical_comp, dft_exch))
            iqa_results['x_hf_total'], iqa_results['coul_total'], \
                iqa_results['x_hf_atomic'], iqa_results['coul_atomic'] = self.ee_iqa_hf()
            iqa_results['xc_total'], iqa_results['xc_atomic'] = self.ee_iqa_dft(rho, dft_xc_edens["edens_xc"])
            if 'coeff_mix' in dft_xc_edens.keys():
                if dft_xc_edens['coeff_mix']:
                    # Scaling HF exchange
                    iqa_results['x_hf_total'] = iqa_results['x_hf_total'] * dft_xc_edens[
                        'coeff_mix']
                    iqa_results['x_hf_atomic'] = iqa_results['x_hf_atomic'] * dft_xc_edens[
                        'coeff_mix']
            else:
                iqa_results.pop('x_hf_total')
                iqa_results.pop('x_hf_atomic')
        elif molecule.lot == 'rhf':
            iqa_results['x_total'], iqa_results['coul_total'], iqa_results['x_atomic'], iqa_results[
                'coul_atomic'] = self.ee_iqa_hf()
        else:
            raise ValueError(f'Need to specify an exchange functional too. Got {dft_exch}')
        logging.info('Summary IQA')
        print('Nucleus-Nucleus repulsion energy: ', iqa_results['nn_total'])
        print('Nucleus-Nucleus repulsion energy (HORTON): ', nn_horton)
        diff = abs(iqa_results['nn_total'] - nn_horton)
        print('DIFF (kcal/mol) ', diff / 0.0015936014376406278)

        print('------------------------------------')
        print('Electron-Nucleus attraction energy: ', iqa_results['en_total'])
        print('Electron-Nucleus attraction energy (HORTON): ', ne_horton)
        diff = abs(iqa_results['en_total'] - ne_horton)
        print('DIFF (kcal/mol) ', diff / 0.0015936014376406278)
        print()
        if part:
            print('Atomic Electron-Nucleus attraction energy:')
            for idx, at in enumerate(molecule.atnums):
                print(f"{at}   {iqa_results['en_atomic'][idx]}")
        print('------------------------------------')
        print('Kinetic energy:  ', iqa_results['kin_total'])
        print('Kinetic energy(HORTON):  ', kinetic_horton)
        diff = abs(iqa_results['kin_total'] - kinetic_horton)
        print('DIFF (kcal/mol) ', diff / 0.0015936014376406278)
        print()
        if part:
            print('Atomic Kinetic energy:')
            for idx, at in enumerate(molecule.atnums):
                print(f"{at}   {iqa_results['kin_atomic'][idx]}")
        print('------------------------------------')
        print('Electron-Electron interaction energy')
        print()
        print('Coulomb energy: ', iqa_results['coul_total'])
        print()
        if part:
            print('Atomic Coulomb energy:')
            for idx, at in enumerate(molecule.atnums):
                print(f"{at}   {iqa_results['coul_atomic'][idx]}")
        print()
        if dft_exch and dft_corr:
            print('DTF Exchange energy: ', iqa_results['x_total'])
            print('DTF Correlation energy: ', iqa_results['c_total'])
            print()
            if part:
                print('Atomic DFT Exchange energy:')
                for idx, at in enumerate(molecule.atnums):
                    print(f"{at}   {iqa_results['x_atomic'][idx]}")
                print('Atomic DFT Correlation energy:')
                for idx, at in enumerate(molecule.atnums):
                    print(f"{at}   {iqa_results['c_atomic'][idx]}")

        if dft_exch and not dft_corr:
            print('DTF Exchange-Correlation energy: ', iqa_results['xc_total'])
            print()
            if part:
                print('Atomic DFT Exchange-Correlation energy:')
                for idx, at in enumerate(molecule.atnums):
                    print(f"{at}   {iqa_results['xc_atomic'][idx]}")

        if 'x_hf_total' in iqa_results.keys():
            print('HF Exchange: ', iqa_results['x_hf_total'])
            print('COEFF MIX :', dft_xc_edens['coeff_mix'])
            print()
            if part:
                print('Atomic HF Exchange energy:')
                for idx, at in enumerate(molecule.atnums):
                    print(f"{at}   {iqa_results['x_hf_atomic'][idx]}")

        total_sum_energy = 0
        for k in iqa_results.keys():
            if k.endswith('total'):
                total_sum_energy += iqa_results[k]

        print('TOTAL ENERGY (SUM OF COMPONENTS):', total_sum_energy)
        print('TOTAL ENERGY (FCHK):', molecule.energy)
        print('DIFF (AU) :', total_sum_energy - molecule.energy)
        print('DIFF (Kcal) :', (total_sum_energy - molecule.energy) / 0.0015936014376406278)

        return iqa_results

    def iqa_pairwise(self, dft_exch=None, dft_corr=None):
        r"""Compute pairwise interaction for exchange-correlation energy density from DFT
        functionals using BOD partition method.
        """

        molecule = self.molecule
        grid = self.grid

        analytical_comp = get_horton_analytical_components(molecule)
        dft_xc_edens = {}
        # assuming dft_corr and dft_exch specified together
        if dft_corr and dft_exch:
            data_meanfield = get_libxc_xc_density(molecule, grid, analytical_comp, dft_exch, dft_corr)
            meanfield_obj_x, meanfield_obj_c  = data_meanfield[-2:]
            dft_exch_ab = self.ee_iqa_dft_pairwise(meanfield_obj_x)
            dft_coul_ab = self.ee_iqa_dft_pairwise(meanfield_obj_c)
            logging.info("CALCULATION: DFT EXCHANGE/CORRELATION PAIRWISE")
            print('dft_exch_ab')
            print(dft_exch_ab)
            print('dft_coul_ab')
            print(dft_coul_ab)
            print('dft_total_ab')
            print(dft_exch_ab + dft_coul_ab)
            logging.info("CALCULATION: HF COULOMB AND EXCHANGE PAIRWISE")
            logging.warning("6n integrals can be long")
            ab_hf_coul, ab_hf_exch = self.ee_iqa_hf_pairwise()
            print('ab_hf_coul')
            print(ab_hf_coul)
            print('ab_hf_exch')
            print(ab_hf_exch)
            iqa_pairwise_results = {
                # 'dft_exch_ab':dft_exch_ab, 'dft_coul_ab':dft_coul_ab,
                 'ab_hf_exch':ab_hf_exch, 'ab_hf_coul':ab_hf_coul
            }

            return iqa_pairwise_results


        elif dft_exch:
            data_meanfield = get_libxc_xc_density(molecule, grid, analytical_comp, dft_exch)
            meanfield_obj_xc = data_meanfield[-1]
            dft_xc_ab = self.ee_iqa_dft_pairwise(meanfield_obj_xc)
            print('dft_xc_ab')
            print(dft_xc_ab)
            ab_hf_coul, ab_hf_exch = self.ee_iqa_hf_pairwise()
            print('ab_hf_coul')
            print(ab_hf_coul)
            print('ab_hf_exch')
            print(ab_hf_exch)
            iqa_pairwise_results = {
                # 'dft_xc_ab': dft_xc_ab,
                'ab_hf_exch': ab_hf_exch, 'ab_hf_coul': ab_hf_coul
            }

            return iqa_pairwise_results



    def nn_iqa(self):
        r"""Compute nuclear-nuclear repulsion energy.

        math::
            \sum_{A>B} \frac{Z_{A}Z_{B}}{R_{A}-R_{B}}

        Where Z_{A} and Z_{B} are nuclear charges.

        Returns
        -------
        total_nn: np.array()
            Total value for nuclear-nuclear repulsion energy
        """

        molecule = self.molecule

        # Compute Nucleus-Nucleus repulsion
        logging.info("CALCULATING NUCLEUS-NUCLEUS REPULSION ENERGY")
        rab = np.triu(np.linalg.norm(molecule.atcoords[:, None] - molecule.atcoords, axis=-1))
        atomic_charges = np.triu(molecule.atnums[:, None] * molecule.atnums)[np.where(rab > 0)]
        rab = rab[rab > 0]
        total_nn = np.sum(atomic_charges / rab)
        print("TOTAL NN ENERGY: ", total_nn)
        print()

        return total_nn

    def en_iqa(self, share_factor=0.5):
        r"""Compute IQA's electron-nuclear attraction energy.

         math::
            \sum_{A, B}\int_{A} \rho(r_{1}) \frac{Z_{B}}{r_{1}-R_{B}} dr_{1}

        Where \rho(r) is the electron density and Z_{B} are nuclear charges.

        Parameters
        ----------
        share_factor: float
            How much of the energy is given to atom i when its density interacts with other nuclei j(x)
            and how much energy is given to atom i when its nucleus interacts with other atomic densities j(1-x)

        Returns
        -------
        total_en: np.array()
            Total value for electron-nuclear attraction energy
        en_cond_en: np.array(natoms)
            Atomic condensed electron-nuclear attraction energies

        """

        molecule = self.molecule
        grid = self.grid
        part = self.part
        dens = self.dens

        # Compute Electron-Nucleus attraction
        logging.info("CALCULATING ELECTRON-NUCLEI ATTRACTION ENERGY")
        natoms = molecule.atnums.shape[0]

        rij = np.linalg.norm(molecule.atcoords[:, None, :] - grid.points, axis=-1)
        total_en = grid.integrate(np.sum((-molecule.atnums[:, None] * (dens / rij)), axis=0))
        print("TOTAL EN ENERGY: ", total_en)

        en_cond_en = None
        if part:
            if part.__class__.__name__ in ['VarHirshfeld', 'HirshfeldI','Hirshfeld']:
                at_weights = part.weights
            else:
                # Doing a for loop to get all at_weights. Using Part object from Horton does not allow
                # to get all at the same time
                at_weights = np.zeros((natoms, grid.points.shape[0]))
                start = 0
                stop = 0
                for i in range(molecule.natom):
                    if part.part.local:
                        stop +=  part.part.cache.load("at_weights", i).shape[0]
                        at_weights[i,start:stop] = part.part.cache.load("at_weights", i)
                        start = stop
                    else:
                        at_weights[i] = part.part.cache.load("at_weights", i)
            # math
            en_atomic = -molecule.atnums[None, :, None] * (
                (dens[None:,] * at_weights)[:, None, :] / rij[None, :, :]
            )
            en_atomic_matrix = np.zeros((natoms, natoms))
            for i in range(natoms):
                for j in range(natoms):
                    if not part.__class__.__name__ in ['VarHirshfeld', 'HirshfeldI', 'Hirshfeld']:
                        if part.part.local:
                            at_grid = part.part.get_grid(i)
                            local_prop = part.part.to_atomic_grid(i, en_atomic[i][j])
                            en_atomic_matrix[i][j] = at_grid.integrate(local_prop)
                        else:
                            en_atomic_matrix[i][j] = grid.integrate(en_atomic[i][j])
                    else:
                        en_atomic_matrix[i][j] = grid.integrate(en_atomic[i][j])

            logging.info("Decomposing Electron-Nuclei energy into atomic contributions")
            print(en_atomic_matrix)
            # share_matrix
            share_factor = share_factor
            logging.info(
                "Condensed into purely Atomic Contributions with "
                f"Sharing Factor x = {share_factor} "
            )
            share_matrix = np.zeros((natoms, natoms, natoms))
            for i in range(natoms):
                share_matrix[i, i, :] = share_factor
                share_matrix[i, :, i] = 1 - share_factor
                share_matrix[i, i, i] = 1

            en_cond_en = np.sum(en_atomic_matrix[None, :, :] * share_matrix, axis=(1, 2))
            print(en_cond_en)
            print(np.sum(en_cond_en), total_en)
            # assert 5 == 6

            #assert_almost_equal(np.sum(en_cond_en), total_en, decimal=2)

        print()
        return total_en, en_cond_en

    def kin_iqa(self):
        r"""Compute IQA's kinetic energy.

        math::
            T_+ (\mathbf{r}_n) = \frac{1}{2} \left. \nabla_{\mathbf{r}}^{2} \gamma(\mathbf{r}, \mathbf{r}')
                                 \right|_{\mathbf{r} = \mathbf{r}' = \mathbf{r}_n}
            Kinetic= T_{\alpha} (\mathbf{r}_n) = T_+(\mathbf{r}_n) + \alpha \nabla^2 \rho(\mathbf{r}_n)

        Where $T_+$ is the positive definite kinetic energy density, with $\gamma$ being the one-electron
        density matrix, and Kinetic is the expression use to return `total_en`.

        Returns
        -------
        total_kin: np.array()
            Total value for kinetic energy
        at_kin: np.array(natoms)
            Atomic kinetic energies.

        """

        molecule = self.molecule
        grid = self.grid
        part = self.part
        basis = self.basis
        dm = self.dm

        logging.info("CALCULATING KINETIC ENERGY")

        natoms = molecule.atnums.shape[0]
        output = evaluate_general_kinetic_energy_density(dm, basis, grid.points, -0.25)
        output_posdef = evaluate_posdef_kinetic_energy_density(dm, basis, grid.points)
        total_kin = grid.integrate(output)
        total_kin_posdef = grid.integrate(output_posdef)
        print("TOTAL KINETIC ENERGY: ", total_kin)

        at_kin = None
        if part:
            if part.__class__.__name__ in ['VarHirshfeld', 'HirshfeldI', 'Hirshfeld']:
                at_weights = part.weights
                at_kin_raw = at_weights * output[None, :]
                at_kin_raw_posdef = at_weights * output_posdef[None, :]
                at_kin = np.array([grid.integrate(at_kin_raw[i]) for i in range(natoms)])
                at_kin_posdef = np.array([grid.integrate(at_kin_raw_posdef[i]) for i in range(natoms)])
                # at_weights = part.weights
                # at_kin = np.array(natoms)
                # for i, at_w in enumerate(part.weights):
                #     at_kin_raw_at = at_w * output
                #     at_kin[i] = grid.integrate(at_kin_raw_at)
            else:
                at_kin = part.condense_to_atoms(output)
                at_kin_posdef = part.condense_to_atoms(output_posdef)
            logging.info("Decomposing Kinetic energy into atomic contributions.")
            print(at_kin)
            # assert 5 == 6

        print()
        return total_kin, at_kin, total_kin_posdef, at_kin_posdef

    def ee_iqa_hf(self):
        r"""Compute Hartree Fock electron-electron interaction energy.

        .. math::
        E_c = \int\int \frac{\varphi^{*}_{i}(r_{1})\varphi^{*}_{j}\varphi_{i}(r_{1})\varphi_{j}(r_{2}) }{r_{1}-r_{2}} dr_{1} dr_{2} = \\
              \int\int \frac{\rho(r_{1}) \rho(r_{2})}{r_{1}-r_{2}} dr_{1} dr_{2} = \\
              \int \rho(r_{1}) dr_{1}\int \frac{\rho(r_{2})}{r_{1}-r_{2}} dr_{2} = \\
              \int \rho(r_{1}) dr_{1} \sum_{\mu\nu} P_{\mu\nu} Vab_{\mu\nu,r_{2}} dr_2 \\
              Vab_{\mu\nu,r_{2}} = \int \frac{\phi^{*}_{\mu}(r_{2})\phi_{\nu}(r_{2})}{r_{1}-r_{2}}

        E_exc =  \int\int \frac{\varphi^{*}_{i}(r_{1})\varphi^{*}_{j}(r_{2})\varphi_{i}(r_{2})\varphi_{j}(r_{1}) }{r_{1}-r_{2}} dr_{1} dr_{2} = \\
                 \int \varphi^{*}_{i}(r_{1})\varphi_{j}(r_{1}) dr_{1} \int  \frac{\varphi^{*}_{j}(r_{2})\varphi_{i}(r_{2})}{r_{1}-r_{2}} dr_{2} \\
                  \int  \frac{\varphi^{*}_{j}(r_{2})\varphi_{i}(r_{2})}{r_{1}-r_{2}} dr_{2} = \sum_{\mu\nu} C\tran_{\mu i} Vab_{\mu\nu,r_{2}} C_{\nu j}

        Where:
            \varphi_{j}(r_{1}): Molecular orbital
            \phi_{\nu}(r_{1}):  Atomic orbital
            Vab_{\mu\nu,r_{2}}: Point charge integrals in the basis of Atomic orbitals with 1 point charge values
                                in the grid points.

        Parameters
        ----------
        dens: np.array(npoints)
            Electron density evaluated at the grid points.

        Returns
        -------
        total_exch: np.array()
            Total value for hartree fock exchange energy
        total_col: np.array()
            Total value for Coulomb energy
        at_exch: np.array(natoms)
            Atomic exchange energies.
        at_colomb: np.array(atoms)
            Atomic exchange energies.
        """

        molecule = self.molecule
        grid = self.grid
        part = self.part
        basis = self.basis
        dm = self.dm
        dens = self.dens

        logging.info("CALCULATING COULOMB AND HF EXCHANGE ENERGY")

        natoms = molecule.atnums.shape[0]
        nao = dm.shape[0]
        eval_ao = evaluate_basis(basis, grid.points)
        # The broadcast operation is the same as the for loop. Summing over mu because broadcasting
        # allocates too much memory
        eval_mo = np.zeros((molecule._iodata.mo.coeffs.T.shape[0], eval_ao.shape[1]))
        for i in range(molecule._iodata.mo.coeffs.T.shape[0]):
            mo = np.zeros((eval_ao.shape[1]))
            for mu in range(eval_ao.shape[0]):
                mo += molecule._iodata.mo.coeffs.T[i, mu] * eval_ao[mu, :]
            eval_mo[i] = mo
        logging.warning("Calculating Coulomb and Exchange: expect long integrals")
        # math
        total_coul_raw = np.zeros(grid.npoints)
        total_exch_raw = np.zeros(grid.npoints)
        istart = 0
        chunk_size = 250000000 // (nao**2)
        while istart < grid.npoints:
            iend = min(istart + chunk_size, grid.npoints)
            # Getting corresponding sliced electron density
            logging.info(
                "Computing for Gridpoint Chunk [ CHUNK START ... END / TOTAL GRIDPOINTS ]: {} ... {} / {}".format(
                    istart, iend, grid.npoints
                )
            )
            chunk_grid = grid.points[istart:iend, :]
            chunk_dens = dens[istart:iend]
            vab = point_charge_integral(basis, chunk_grid, np.ones(chunk_dens.shape[0]))
            # Coulomb
            vab_rho = np.trace(np.tensordot(dm, (vab * chunk_dens), axes=(1, 0)))
            total_coul_chunk = -0.5 * vab_rho
            # Exchange
            # Because only restricted _occs_a.shape[0] == _occs_b.shape[0]
            occupied_mo = np.zeros(molecule._mo._occs_a.shape[0])
            occupied_mo[molecule._mo._occs_a> 0] = 1
            t1 = np.einsum("abn,ai->ibn", vab, molecule._iodata.mo.coeffs)
            t1 = np.einsum("ibn,bj->ijn", t1, molecule._iodata.mo.coeffs)
            t1 = t1 * occupied_mo[:, None, None]
            total_exch_chunk = np.einsum(
                "ijn,in,jn->n",
                t1,
                (eval_mo[:, istart:iend] * occupied_mo[:, None]),
                (eval_mo[:, istart:iend] * occupied_mo[:, None]),
            )
            del t1

            total_coul_raw[istart:iend] = total_coul_chunk
            total_exch_raw[istart:iend] = total_exch_chunk
            istart = iend

        total_coul = grid.integrate(total_coul_raw)
        total_exch = grid.integrate(total_exch_raw)

        print("TOTAL COULOMB: ", total_coul)
        print("TOTAL EXCHANGE: ", total_exch)

        at_exch = None
        at_coulomb = None
        if part:
            if part.__class__.__name__ in ['VarHirshfeld', 'HirshfeldI','Hirshfeld']:
                at_weights = part.weights
                at_coulomb_raw = at_weights * total_coul_raw[None, :]
                at_exch_raw = at_weights * total_exch_raw[None, :]
                at_coulomb = np.array([grid.integrate(at_coulomb_raw[i]) for i in range(natoms)])
                at_exch = np.array([grid.integrate(at_exch_raw[i]) for i in range(natoms)])
            else:
                at_coulomb = part.condense_to_atoms(total_coul_raw)
                at_exch = part.condense_to_atoms(total_exch_raw)
            logging.info("Decomposing Coulomb and Exchange into atomic contributions.")
            print("Coulomb")
            print(at_coulomb, np.sum(at_coulomb))
            print("Exchange")
            print(at_exch, np.sum(at_exch))

        print()
        return total_exch, total_coul, at_exch, at_coulomb

    def ee_iqa_hf_pairwise(self):
        r"""Compute Hartree Fock electron-electron interaction energy pairwise interactions.


        Returns
        -------
        total_exch_aa: np.array()
            AA value for hartree fock exchange energy
        total_exch_ab: np.array()
            AB value for hartree fock exchange energy
        total_col_aa: np.array()
            AA value for Coulomb energy
        total_col_aa: np.array()
            AB value for Coulomb energy


        Note: Sum of AA and AB terms should get back ee_iqa_hf
        """

        # if not isinstance(self.grid._grid, MolGrid) or not isinstance(self.grid_2._grid, MolGrid):
        #     raise TypeError(
        #         f"Both grids should be an instance of MolGrid from qc-devs."
        #     )

        # Draft 6N integration
        molecule = self.molecule
        basis = self.basis
        part = self.part
        part2 = self.part_2
        grid1 = self.grid
        grid2 = self.grid_2
        points1 = self.grid._grid.points
        points2 = self.grid_2._grid.points
        # Grid1
        eval_ao_g1 = evaluate_basis(basis, points1)
        eval_mo_g1 = compute_molecular_orbitals_from_ao(molecule, eval_ao_g1)
        # Grid2
        eval_ao_g2 = evaluate_basis(basis, points2)
        eval_mo_g2 = compute_molecular_orbitals_from_ao(molecule, eval_ao_g2)

        # Compute density/exchange density from Molecular orbitals
        # Because only restricted _occs_a.shape[0] == _occs_b.shape[0]
        occupied_mo = molecule._iodata.mo.occs[molecule._iodata.mo.occs > 0].shape[0]

        if part:
            if part2:
                logging.warning("Molecular 6N integration using local grids")
                ab_hf_coul = np.zeros((len(self.molecule.numbers), len(self.molecule.numbers)))
                ab_hf_exch = np.zeros((len(self.molecule.numbers), len(self.molecule.numbers)))
                natoms = len(self.molecule.numbers)
            else:
                raise ValueError("6N integration uses local grids and 2 partition objects")

        for a in range(natoms):
            if isinstance(self.grid._grid, MolGrid):
                p_start_a = grid1.indices[a]
                p_end_a = grid1.indices[a + 1]
                atgrid_a = self.grid._grid._atgrids[a]
            elif self.grid._grid.__class__.__name__ == 'BeckeMolGrid':
                atgrid_a = grid1.subgrids[a]
                p_start_a = atgrid_a.begin
                p_end_a = atgrid_a.end
            for b in range(a + 1):
                if isinstance(self.grid_2._grid, MolGrid):
                    p_start_b = grid2.indices[b]
                    p_end_b = grid2.indices[b + 1]
                    atgrid_b = grid2._atgrids[b]
                elif self.grid_2._grid.__class__.__name__ == 'BeckeMolGrid':
                    atgrid_b = grid2.subgrids[b]
                    p_start_b = atgrid_b.begin
                    p_end_b = atgrid_b.end
                integral_coul_ab = 0
                integral_ex_ab = 0
                for i_d1 in range(atgrid_a.points.shape[0]):
                    progress = i_d1 / atgrid_a.points.shape[0]
                    if progress * 100 in [15.0, 25.0, 50.0, 75.0, 90.0]:
                        # print('asdfasdf')
                        print('Progress: ', progress * 100)
                    d_coul = 0
                    d_ex = 0
                    for i in range(occupied_mo):
                        for j in range(occupied_mo):
                            d_coul += 2 * eval_mo_g1[i, i_d1 + p_start_a] * eval_mo_g1[i, i_d1 + p_start_a] * eval_mo_g2[j] * eval_mo_g2[j]
                            d_ex += eval_mo_g1[i, i_d1 + p_start_a] * eval_mo_g1[j, i_d1 + p_start_a] * eval_mo_g2[j] * eval_mo_g2[i]

                    rij = np.linalg.norm(atgrid_a.points[i_d1, :] - points2, axis=-1)
                    rij[rij == 0] = 1.0e-9
                    d_ex = d_ex / rij
                    d_coul = d_coul / rij
                    # Subset part weights b
                    w_subset_b = part2.weights[b][p_start_b: p_end_b]
                    # Subset d_coul and d_ex
                    d_coul_b = d_coul[p_start_b: p_end_b]
                    d_ex_b = d_ex[p_start_b: p_end_b]
                    # Atomic integration for r2 using part weights
                    part_coul_b = atgrid_b.integrate(d_coul_b * w_subset_b)
                    part_ex_b = atgrid_b.integrate(d_ex_b * w_subset_b)
                    integral_coul_ab += part_coul_b * (atgrid_a.weights[i_d1] * part.weights[a, p_start_a : p_end_a][i_d1])
                    integral_ex_ab += part_ex_b * (atgrid_a.weights[i_d1] * part.weights[a, p_start_a : p_end_a][i_d1])

                ab_hf_coul[a, b] = integral_coul_ab
                ab_hf_exch[a, b] = integral_ex_ab
                if a != b:
                    ab_hf_coul[b, a] = integral_coul_ab
                    ab_hf_exch[b, a] = integral_ex_ab


        return ab_hf_coul, ab_hf_exch

    def ee_iqa_dft(self, dft_dens, dft_xc_dens):

        grid = self.grid
        part = self.part
        molecule = self.molecule

        natoms = molecule.atnums.shape[0]

        dft_xc_total = grid.integrate(dft_xc_dens, dft_dens)

        at_dft_xc = None
        if part:
            if part.__class__.__name__ in ['VarHirshfeld', 'HirshfeldI','Hirshfeld']:
                at_weights = part.weights
                at_dft_xc_raw = at_weights * dft_xc_dens[None, :]
                at_dft_xc = np.array([grid.integrate(at_dft_xc_raw[i], dft_dens) for i in range(natoms)])
            else:
                at_dft_xc = part.condense_to_atoms((dft_xc_dens * dft_dens))
            logging.info("Decomposing XC into atomic contributions.")
            print('XC component')
            print(at_dft_xc)
            print(np.sum(at_dft_xc), dft_xc_total)
            assert_almost_equal(np.sum(at_dft_xc), dft_xc_total, decimal=2)

        print()
        return dft_xc_total, at_dft_xc
    # def ee_iqa_dft_pairwise(self, dft_dens, dft_xc_dens):
    def ee_iqa_dft_pairwise(self, meanfield_dft):

        grid = self.grid
        part = self.part
        molecule = self.molecule
        basis = self.basis

        # Get different data
        at_weights = part.weights
        natoms = molecule.atnums.shape[0]


        # BOD for LDA
        # Evaluate Atomic orbitals basis
        eval_ao = evaluate_basis(basis, grid.points)
        # Get Molecular orbitals from AO
        eval_mo = compute_molecular_orbitals_from_ao(molecule, eval_ao)
        # Computing Molecular overlap from eval_mo
        sij = eval_mo[None, :, :] * eval_mo[:, None, :]
        # Getting atomic overlap matrices
        sij_at_1 = at_weights[:, None, None, :] * sij[None, :, :, :]

        # Computing BOD for each pair of atoms
        # print(natoms, natoms, grid.points.shape[0])
        bod = np.zeros((natoms, natoms, grid.points.shape[0]))
        for a in range(natoms):
            for b in range(natoms):
                for i in range(molecule._iodata.mo.occs[molecule._iodata.mo.occs > 0].shape[0]):
                    for j in range(molecule._iodata.mo.occs[molecule._iodata.mo.occs > 0].shape[0]):
                            bod[a,b] += 2* (at_weights[a] * grid.integrate(sij_at_1[b][i][j]) +
                                            at_weights[b] * grid.integrate(sij_at_1[a][j][i])) * \
                                            eval_mo[i] * eval_mo[j]
        # Compute Grad BOD for GGA and Hyb
        if meanfield_dft.__class__.__name__ in ["RLibXCGGA", "RLibXCHybridGGA"]:
            # BOD for GGA and Hybrid GGA
            # Evaluate molecular orbitals gradient
            # todo: this could potentially go to utils but then Gbasis needs to be imported there
            orders_d = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            eval_mo_grad = []
            eval_mo_grad_condensed = []
            for ind, orders in enumerate(orders_d):
                eval_mo_grad_order = evaluate_deriv_basis(basis, grid.points, np.array(orders),
                                                             transform=molecule._iodata.mo.coeffs.T)
                # print(eval_mo_grad_order.shape)
                eval_mo_grad.append(eval_mo_grad_order[:, :, None])
                eval_mo_grad_condensed.append(eval_mo_grad_order)
            eval_mo_grad = np.concatenate(eval_mo_grad, axis=2)

            # Getting the norm-squared of the gradient of the spin-summed electron density(sigma_full)
            # grad = evaluate_density_gradient(self.dm, basis, grid.points)
            # sigma_full = 4*((grad/2)**2).sum(axis=1)

            # print(sigma_full)

            # Computing BOD for GGA functionals for each pair of atoms
            print(natoms, natoms, grid.points.shape[0])
            grad_bod = np.zeros((natoms, natoms, grid.points.shape[0]))
            for a in range(natoms):
                for b in range(natoms):
                    for i in range(molecule._iodata.mo.occs[molecule._iodata.mo.occs > 0].shape[0]):
                        for j in range(molecule._iodata.mo.occs[molecule._iodata.mo.occs > 0].shape[0]):
                            grad_bod[a, b] += (np.sum((at_weights[a, :, None]* grid.integrate(sij_at_1[b][i][j]) + at_weights[b, :, None]* grid.integrate(sij_at_1[a][i][j])) * eval_mo_grad[i, :, :] * eval_mo[j, :,  None], axis=-1))


            grad_bod = 16 * (grad_bod ** 2)

        # Integration BOD
        int_bod = 0
        for a in range(natoms):
            for b in range(a):
                int_bod += grid.integrate(bod[a,b])
        print('Integration BOD (TOTAL)', int_bod)


        # Combine atomic partition weights
        combine_atweights = np.zeros((at_weights[0].shape[0]))
        # Below code for grid is for Horton grid.
        # todo: add code for qc-devs grid
        for at in range(natoms):
            start = grid._grid.subgrids[at].begin
            end =  grid._grid.subgrids[at].end
            combine_atweights[start:end] = at_weights[at][start:end]

        # LDA
        if meanfield_dft.__class__.__name__ == "RLibXCLDA":
            dft_pairwise = np.zeros((natoms, natoms))
            for a in range(natoms):
                for b in range(natoms):
                    cache = Cache()
                    cache['rho_full'] = bod[a,b]
                    cache['combine_weights'] = combine_atweights
                    dft_pairwise[a,b] = meanfield_dft.compute_energy(cache, grid._grid)

        #GGA and Hyb
        # print(meanfield_dft.__class__.__name__)
        if meanfield_dft.__class__.__name__ in ["RLibXCGGA", "RLibXCHybridGGA"]:
            dft_pairwise = np.zeros((natoms, natoms))
            for a in range(natoms):
                for b in range(natoms):
                    cache = Cache()
                    cache['rho_full'] = bod[a,b]
                    cache['sigma_full'] = grad_bod[a,b]
                    cache['combine_weights'] = combine_atweights
                    dft_pairwise[a, b] = meanfield_dft.compute_energy(cache, grid._grid)


        return dft_pairwise



    # def _coul_pairwise(self, grid1, grid2):
    #     grid1 = grid1.reshape(1, -1)
    #     rij = np.linalg.norm(grid1 - grid2, axis=-1)
    #     density_1 = evaluate_density(self.dm, self.basis, grid1)
    #     density_2 = evaluate_density(self.dm, self.basis, grid2)
    #     rij[rij == 0] = 1.0e-9
    #     output = density_1 * (density_2 / rij)
    #     return output
    #
    # def _exch_pairwise(self, grid1, grid2):
    #
    #     molecule = self.molecule
    #     basis = self.basis
    #     part = self.part
    #     # self.part.weights[0], self.part_2.weights[1]
    #     grid1 = grid1.reshape(1, -1)
    #     # Evaluate basis set orbitals
    #     # grid1 == 1 point
    #     eval_ao_g1 = evaluate_basis(basis, grid1)
    #     # grid2 == Npoints
    #     eval_ao_g2 = evaluate_basis(basis, grid2)
    #     # Convert to Molecular orbitals
    #     eval_mo_g1 = np.zeros((molecule._iodata.mo.coeffs.T.shape[0], eval_ao_g1.shape[1]))
    #     eval_mo_g2 = np.zeros((molecule._iodata.mo.coeffs.T.shape[0], eval_ao_g2.shape[1]))
    #     for i in range(molecule._iodata.mo.coeffs.T.shape[0]):
    #         mo_g1 = np.zeros((eval_ao_g1.shape[1]))
    #         mo_g2 = np.zeros((eval_ao_g2.shape[1]))
    #         for mu in range(eval_ao_g1.shape[0]):
    #             mo_g1 += molecule._iodata.mo.coeffs.T[i, mu] * eval_ao_g1[mu, :]
    #             mo_g2 += molecule._iodata.mo.coeffs.T[i, mu] * eval_ao_g2[mu, :]
    #         eval_mo_g1[i] = mo_g1
    #         eval_mo_g2[i] = mo_g2
    #
    #     # Compute exchange density from Molecular orbitals
    #     # Because only restricted _occs_a.shape[0] == _occs_b.shape[0]
    #     occupied_mo = molecule._iodata.mo.occs[molecule._iodata.mo.occs > 0].shape[0]
    #     d_ex = 0
    #     for i in range(occupied_mo):
    #         for j in range(occupied_mo):
    #             d_ex += eval_mo_g1[i] * eval_mo_g1[j] * eval_mo_g2[j] * eval_mo_g2[i]
    #             print('d_ex', d_ex)
    #
    #     rij = np.linalg.norm(grid1 - grid2, axis=-1)
    #     rij[rij == 0] = 1.0e-9
    #     output = d_ex / rij
    #
    #     return -output
