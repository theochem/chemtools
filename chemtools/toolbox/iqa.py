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

import numpy as np
from numpy.testing import assert_almost_equal

from chemtools.toolbox.utils import check_arg_molecule, \
                                    get_horton_analytical_components, \
                                    get_molecular_grid, get_libxc_xc_density
from iodata import load_one
from gbasis.wrappers import from_iodata
from gbasis.evals.density import evaluate_density
from gbasis.evals.density import evaluate_general_kinetic_energy_density
from gbasis.evals.eval import evaluate_basis
from gbasis.integrals.point_charge import point_charge_integral
from chemtools.utils.cube import UniformGrid
from chemtools.wrappers.grid import MolecularGrid
from chemtools.wrappers.molecule import Molecule
from chemtools.wrappers.part import DensPart


class IQA(object):
    """Interacting Quantum Atoms (IQA) Class."""

    def __init__(self, molecule, basis, dm, grid, part=None, molecule_chemtools=None):
        """Initialize class.

        Parameters
        ----------
        molecule : `Molecule`
            Instance of `Molecular` class from IOData.
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
        molecule_chemtools: `Molecule`, optional
            Instance of `Molecular` class from Chemtools(Horton). This molecule class is needed
            to perform IQA decomposition for DFT wavefunctions.

        """
        # check basis
        if not basis[0][0].__class__.__name__ == "IODataShell":
            raise TypeError("basis should have been created with Gbasis")

        # check dm
        if not (isinstance(dm, np.ndarray) and dm.ndim == 2 and dm.dtype == float):
            raise TypeError("One-electron density matrix must be a 2D-array float `dtype`")
        if dm.shape[0] != dm.shape[1]:
            raise ValueError("One-electron density matrix must be a square matrix.")
        if not np.allclose(dm, dm.T):
            raise ValueError("One-electron density matrix must be symmetric.")

        # check grid and part
        if part is not None and not isinstance(part, DensPart):
            raise TypeError("Argument part should be an instance of DensPart class.")
        if part is not None and not (part.numbers == molecule.atnums).all():
            raise ValueError("DensPart molecule different from molecule")
        if not isinstance(grid, UniformGrid) and not isinstance(grid, MolecularGrid):
            raise TypeError(
                "Argument part should be an instance of MolecularGrid or UniformGrid class."
            )
        if not (grid.numbers == molecule.atnums).all():
            raise ValueError("Grid molecule different from molecule")

        self.molecule = molecule
        if molecule_chemtools:
            self.molecule_chemtools = molecule_chemtools
        self.basis = basis
        self.dm = dm
        self.grid = grid
        self.part = None
        if part is not None:
            self.part = part

    @classmethod
    def from_molecule(cls, molecule_iodata, molecule_chemtools, grid, part=None, scheme=None):
        """Initialize Interacting Quantum Atoms (IQA) class from `Molecules` object.

        Parameters
        ----------
        molecule_iodata : `Molecule`
            Instance of `Molecular` class from IOData.
        molecule_chemtools : `Molecule`
             Instance of `Molecular` class from Chemtools(Horton).
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
        molecule_chemtools = check_arg_molecule(molecule_chemtools)
        # Check molecule_iodata
        if not (
            molecule_iodata.__class__.__name__ == "IOData" and hasattr(molecule_iodata, "obasis")
        ):
            raise ValueError("`molecule_iodata` must be a `IOData` instance.")
        if not (molecule_chemtools.numbers == molecule_iodata.atnums).all():
            raise ValueError(
                "molecule_chemtools and molecule_iodata do not correspond to the same calculation"
            )
        # Check restricted closed shell single-determinant wave-function (molecule_iodata)
        if molecule_iodata.mo.kind != "restricted":
            raise ValueError("Currently, only 'restricted' wave-functions are supported")
        if int(np.sum(molecule_iodata.mo.occsa)) != int(np.sum(molecule_iodata.mo.occsb)):
            raise ValueError("Code is not tested for open-shell wave-functions.")

        # Check Grid
        if not isinstance(grid, UniformGrid) and not isinstance(grid, MolecularGrid):
            raise TypeError(
                f"Argument grid should be an instance of MolecularGrid or UniformGrid class. Got {grid.__class__.__name__}"
            )
        if not (molecule_iodata.atcoords == grid.centers).all():
            raise ValueError("Molecule and Grid initialized from different molecules")

        # Check/Initialize DensPart object
        if part is None and scheme is None:
            print("No atomic partition scheme provided. No atomic decomposition will be performed.")
        elif part is not None:
            if not isinstance(part, DensPart):
                raise TypeError("Argument part should be an instance of DensPart class.")
        elif scheme == 'H':
            part = DensPart.from_molecule(molecule_chemtools, grid=grid, scheme="h",
                                              local=False)
        elif scheme == 'HI':
            part = DensPart.from_molecule(molecule_chemtools, grid=grid, scheme="hi",
                                              local=False)
        elif scheme is not None and scheme not in ['H', 'HI']:
            raise NotImplementedError(f"Atomic partition {scheme} not yet available")

        # Initialize gbasis
        basis = from_iodata(molecule_iodata)
        one_rdm = molecule_iodata.one_rdms.get("post_scf", molecule_iodata.one_rdms.get("scf"))
        # Check if dm present
        if one_rdm is None:
            # Check if mo present to construct dm
            if molecule_iodata.mo is None:
                raise NotImplementedError("Missing molecule.mo object.Density matrix...")
            else:
                print("Couldn't read Density matrix. Calculating from its components.")
                coeffs, occs = molecule_iodata.mo.coeffs, molecule_iodata.mo.occs
                one_rdm = np.dot(coeffs * occs, coeffs.T)
        # check scf one_rdm is the same for molecule_chemtools and molecule_iodata
        assert_almost_equal(one_rdm, molecule_chemtools.mo.compute_dm(), decimal=1)

        return cls(
            molecule_iodata, basis, one_rdm, grid, part, molecule_chemtools=molecule_chemtools
        )

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

        molecule_iodata = load_one(fname)
        molecule_chemtools = Molecule.from_file(fname)

        # Initialize grid
        # TODO: Add grid type argument to get_molecular_grid
        grid = get_molecular_grid(molecule_chemtools)

        return cls.from_molecule(molecule_iodata, molecule_chemtools, grid, scheme=scheme)

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
        molecule_chemtools = self.molecule_chemtools
        basis = self.basis
        dm = self.dm
        grid = self.grid
        part = self.part

        # Computing electron density
        rho = evaluate_density(dm, basis[0], grid.points, coord_type=basis[1])

        # Initialize results dict
        iqa_results = {}

        logging.info('INITIALIZING INTERACTING QUANTUM ATOMS(IQA) CALCULATION')
        iqa_results['nn_total'] = self.nn_iqa()
        iqa_results['en_total'], iqa_results['en_atomic'] = self.en_iqa(dens=rho)
        iqa_results['kin_total'], iqa_results['kin_atomic'] = self.kin_iqa()
        analytical_comp = get_horton_analytical_components(molecule_chemtools, grid)
        nn_horton = analytical_comp['nn']
        ne_horton = np.trace(dm.dot(analytical_comp['na']))
        kinetic_horton = np.trace(dm.dot(analytical_comp['kin']))
        print('KINETIC DEBUG :', iqa_results['kin_total'])
        print('KINETIC DEBUG (HORTON):', kinetic_horton)
        diff = abs(iqa_results['kin_total'] - kinetic_horton)
        print('DIFF (DEBUG)(kcal/mol) ', diff / 0.0015936014376406278)
        print(grid._grid.points.shape)
        print()
        # todo: test gbasis

        dft_xc_edens = {}
        # assuming dft_corr and dft_exch specified together
        if dft_corr and dft_exch:
            dft_xc_edens.update(get_libxc_xc_density(molecule_chemtools, grid, analytical_comp, dft_exch, dft_corr))
            iqa_results['x_hf_total'], iqa_results['coul_total'], iqa_results['x_hf_atomic'], \
            iqa_results['coul_atomic'] = self.ee_iqa_hf(dens=rho)
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
                get_libxc_xc_density(molecule_chemtools, grid, analytical_comp, dft_exch))
            iqa_results['x_hf_total'], iqa_results['coul_total'], \
                iqa_results['x_hf_atomic'], iqa_results['coul_atomic'] = self.ee_iqa_hf(dens=rho)
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
                'coul_atomic'] = self.ee_iqa_hf(dens=rho)
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

    def en_iqa(self, dens=None, share_factor=0.5):
        r"""Compute IQA's electron-nuclear attraction energy.

         math::
            \sum_{A, B}\int_{A} \rho(r_{1}) \frac{Z_{B}}{r_{1}-R_{B}} dr_{1}

        Where \rho(r) is the electron density and Z_{B} are nuclear charges.

        Parameters
        ----------
        dens: np.array(npoints)
            Electron density evaluated at the grid points.
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
        basis = self.basis

        # Compute Electron-Nucleus attraction
        logging.info("CALCULATING ELECTRON-NUCLEI ATTRACTION ENERGY")
        natoms = molecule.atnums.shape[0]
        if dens is None:
            dm = molecule.one_rdms.get("post_scf", molecule.one_rdms.get("scf"))
            dens = evaluate_density(dm, basis[0], grid.points, coord_type=basis[1])

        rij = np.linalg.norm(molecule.atcoords[:, None, :] - grid.points, axis=-1)
        total_en = grid.integrate(np.sum((-molecule.atnums[:, None] * (dens / rij)), axis=0))
        print("TOTAL EN ENERGY: ", total_en)

        en_cond_en = None
        if part:
            # Doing a for loop to get all at_weights. Using Part object from Horton does not allow
            # to get all at the same time
            at_weights = np.zeros((natoms, grid.points.shape[0]))
            for i in range(molecule.natom):
                at_weights[i] = part.part.cache.load("at_weights", i)
            # math
            en_atomic = -molecule.atnums[None, :, None] * (
                (dens[None:,] * at_weights)[:, None, :] / rij[None, :, :]
            )
            en_atomic_matrix = np.zeros((natoms, natoms))
            for i in range(natoms):
                for j in range(natoms):
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
            # print(np.sum(en_cond_en))

            assert_almost_equal(np.sum(en_cond_en), total_en, decimal=3)

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

        logging.info("CALCULATING KINETIC ENERGY")
        dm = molecule.one_rdms.get("post_scf", molecule.one_rdms.get("scf"))

        natoms = molecule.atnums.shape[0]
        output = evaluate_general_kinetic_energy_density(
            dm, basis[0], grid.points, -0.25, coord_type=basis[1]
        )
        total_kin = grid.integrate(output)
        print("TOTAL KINETIC ENERGY: ", total_kin)

        at_kin = None
        if part:
            # Doing a for loop to get all at_weights. Using Part object from Horton does not allow
            # to get all at the same time
            at_weights = np.zeros((natoms, grid.points.shape[0]))
            for i in range(molecule.natom):
                at_weights[i] = part.part.cache.load("at_weights", i)
            at_kin_raw = at_weights * output[None, :]
            at_kin = np.array([grid.integrate(at_kin_raw[i]) for i in range(natoms)])
            logging.info("Decomposing Kinetic energy into atomic contributions.")
            print(at_kin)

        print()
        return total_kin, at_kin

    def ee_iqa_hf(self, dens=None):
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

        logging.info("CALCULATING COULOMB AND HF EXCHANGE ENERGY")
        dm = molecule.one_rdms.get("post_scf", molecule.one_rdms.get("scf"))
        if dens is None:
            dens = evaluate_density(dm, basis[0], grid.points, coord_type=basis[1])

        natoms = molecule.atnums.shape[0]
        nao = dm.shape[0]
        eval_ao = evaluate_basis(basis[0], grid.points, coord_type=basis[1])
        eval_mo = np.sum((molecule.mo.coeffs.T[:, :, None] * eval_ao), axis=1)
        logging.warning("Calculating Coulomb and Exchange: expect long integrals")
        # Eval falta math
        total_col_raw = np.zeros((grid.npoints))
        total_exch_raw = np.zeros((grid.npoints))
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
            vab = point_charge_integral(
                basis[0], chunk_grid, np.ones(chunk_dens.shape[0]), coord_type=basis[1]
            )
            # Coulomb
            vab_rho = np.trace(np.tensordot(dm, (vab * chunk_dens), axes=(1, 0)))
            total_col_chunk = -0.5 * vab_rho
            # Exchange
            occupied_mo = np.zeros(molecule.mo.occs.shape[0])
            occupied_mo[molecule.mo.occs > 0] = 1
            t1 = np.einsum("abn,ai->ibn", vab, molecule.mo.coeffs)
            t2 = np.einsum("ibn,bj->ijn", t1, molecule.mo.coeffs)
            total_exch_chunk = np.einsum(
                "ijn,in,jn->n",
                (t2 * occupied_mo[:, None, None]),
                (eval_mo[:, istart:iend] * occupied_mo[:, None]),
                (eval_mo[:, istart:iend] * occupied_mo[:, None]),
            )

            total_col_raw[istart:iend] = total_col_chunk
            total_exch_raw[istart:iend] = total_exch_chunk
            istart = iend

        total_col = grid.integrate(total_col_raw)
        total_exch = grid.integrate(total_exch_raw)

        print("TOTAL COULOMB: ", total_col)
        print("TOTAL EXCHANGE: ", total_exch)

        at_exch = None
        at_colomb = None
        if part:
            # Doing a for loop to get all at_weights. Using Part object from Horton does not allow
            # to get all at the same time
            at_weights = np.zeros((natoms, grid.points.shape[0]))
            for i in range(molecule.natom):
                at_weights[i] = part.part.cache.load("at_weights", i)
            at_coulomb_raw = at_weights * total_col_raw[None, :]
            at_exch_raw = at_weights * total_exch_raw[None, :]
            at_colomb = np.array([grid.integrate(at_coulomb_raw[i]) for i in range(natoms)])
            at_exch = np.array([grid.integrate(at_exch_raw[i]) for i in range(natoms)])
            logging.info("Decomposing Coulomb and Exchange into atomic contributions.")
            print("Exchange")
            print(at_colomb)
            print("Coulomb")
            print(at_exch)

        print()
        return total_exch, total_col, at_exch, at_colomb

    def ee_iqa_dft(self, dft_dens, dft_xc_dens):

        molecule = self.molecule
        grid = self.grid
        part = self.part

        natoms = molecule.atnums.shape[0]

        dft_xc_total = grid.integrate(dft_xc_dens, dft_dens)

        at_dft_xc = None
        if part:
            # Doing a for loop to get all at_weights. Using Part object from Horton does not allow
            # to get all at the same time
            at_weights = np.zeros((natoms, grid.points.shape[0]))
            for i in range(molecule.natom):
                at_weights[i] = part.part.cache.load("at_weights", i)
            at_dft_xc_raw = at_weights * dft_xc_dens[None, :]
            at_dft_xc = np.array(
                [grid.integrate(at_dft_xc_raw[i], dft_dens) for i in range(natoms)])
            logging.info("Decomposing XC into atomic contributions.")
            print('XC component')
            print(at_dft_xc)
            assert_almost_equal(np.sum(at_dft_xc), dft_xc_total, decimal=2)

        print()
        return dft_xc_total, at_dft_xc
