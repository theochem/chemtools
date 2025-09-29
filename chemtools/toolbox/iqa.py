# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2024 The ChemTools Development Team
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
from os.path import dirname, join
import os as os
from glob import glob
import multiprocessing as mp
from multiprocessing import set_start_method

import numpy as np
from numpy.testing import assert_almost_equal

import pylibxc as pylibxc

from chemtools.toolbox.utils import (
    check_arg_molecule,
    get_molecular_grid,
    compute_molecular_orbitals_from_ao,
)

import gbasis as gbasis
from gbasis.wrappers import from_iodata
from gbasis.evals.density import evaluate_density
from gbasis.evals.density import evaluate_general_kinetic_energy_density
from gbasis.evals.density import evaluate_posdef_kinetic_energy_density
from gbasis.evals.density import evaluate_density_gradient
from gbasis.evals.eval import evaluate_basis
from gbasis.evals.eval_deriv import evaluate_deriv_basis
from gbasis.integrals.point_charge import point_charge_integral
from gbasis.integrals.kinetic_energy import kinetic_energy_integral
from gbasis.integrals.nuclear_electron_attraction import nuclear_electron_attraction_integral
from gbasis.integrals.electron_repulsion import electron_repulsion_integral

from iodata.periodic import num2sym
from iodata.utils import kcalmol

from grid.basegrid import Grid

from chemtools.utils.cube import UniformGrid
from chemtools.wrappers.grid import MolecularGrid
from chemtools.wrappers.molecule import Molecule
from chemtools.wrappers.part import DensPart


def _integral_point_charge_electron(i, istart, iend, eval_mo, grid, dens, molecule, dm):
    """Auxiliary function for parallelization of coulomb and exchange terms for atomic IQA."""

    from gbasis.wrappers import from_iodata

    basis = from_iodata(molecule)

    logging.info(
        "Computing for Gridpoint Chunk [ CHUNK START ... END / TOTAL GRIDPOINTS ]: {} ... {} / {}".format(
            istart, iend, grid.points.shape[0]
        )
    )
    chunk_grid = grid.points[istart:iend, :]
    chunk_dens = dens[istart:iend]
    vab = point_charge_integral(basis, chunk_grid, np.ones(chunk_dens.shape[0]))
    logging.info("SAFETY CHECK: POINT CHARGE INTEGRAL FINISHED")
    # Coulomb
    vab_rho = np.trace(np.tensordot(dm, (vab * chunk_dens), axes=(1, 0)))
    total_coul_chunk = -0.5 * vab_rho
    # Exchange
    # Because only restricted _occs_a.shape[0] == _occs_b.shape[0]
    occupied_mo = np.zeros(molecule.mo.occsa.shape[0])
    occupied_mo[molecule.mo.occsa > 0] = 1
    t1 = np.einsum("abn,ai->ibn", vab, molecule.mo.coeffs)
    t1 = np.einsum("ibn,bj->ijn", t1, molecule.mo.coeffs)
    t1 = t1 * occupied_mo[:, None, None]
    total_exch_chunk = np.einsum(
        "ijn,in,jn->n",
        t1,
        (eval_mo[:, istart:iend] * occupied_mo[:, None]),
        (eval_mo[:, istart:iend] * occupied_mo[:, None]),
    )
    del t1

    return tuple([istart, iend, total_coul_chunk, total_exch_chunk])


class IQA(object):
    """Interacting Quantum Atoms (IQA) Class."""

    def __init__(
        self, molecule, basis, dm, grid, part=None, grid_2=None, part_2=None, threshold=10
    ):
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

        # check part
        if part is not None:
            if not isinstance(part, DensPart) and part.__class__.__name__ not in [
                "VarHirshfeld",
                "HirshfeldI",
                "Hirshfeld",
            ]:
                raise TypeError(
                    "Argument part should be an instance of DensPart class or VarHirshfeld from rhopart."
                )
        if part is not None and not (part.numbers == molecule.numbers).all():
            raise ValueError("DensPart molecule different from molecule")
        # Check Grid
        if (
            not isinstance(grid, UniformGrid)
            and not isinstance(grid, MolecularGrid)
            and not isinstance(grid, Grid)
        ):
            raise TypeError(
                f"Argument grid should be an instance of MolecularGrid or UniformGrid class. Got {grid.__class__.__name__}"
            )
        # Check optional grid_2 and part_2
        if grid_2:
            if not isinstance(grid_2, UniformGrid) and not isinstance(grid_2, MolecularGrid):
                raise TypeError(
                    "Argument part should be an instance of MolecularGrid or UniformGrid class."
                )
            if not (molecule.coordinates == grid_2.coordinates).all():
                raise ValueError("Molecule and Grid-2 initialized from different molecules")
            if not (grid_2.numbers == molecule.numbers).all():
                raise ValueError("Grid-2 molecule different from molecule")
        if part_2 is not None:
            if not isinstance(part, DensPart) and part.__class__.__name__ not in [
                "VarHirshfeld",
                "HirshfeldI",
                "Hirshfeld",
            ]:
                raise TypeError(
                    "Argument part should be an instance of DensPart class or VarHirshfeld from rhopart."
                )
        if part_2 is not None and not (part.numbers == molecule.numbers).all():
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
        self.integrands, self.totals = self.compare_numerical_analytical(threshold=threshold)

    @classmethod
    def from_molecule(cls, molecule, grid, part=None, scheme=None, grid_2=None, part_2=None):
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
        if not (molecule.coordinates == grid.coordinates).all():
            raise ValueError("Molecule and Grid initialized from different molecules")

        # Check/Initialize DensPart object
        if part is None and scheme is None:
            print("No atomic partition scheme provided. No atomic decomposition will be performed.")
        elif part is not None:
            if not isinstance(part, DensPart):
                raise TypeError("Argument part should be an instance of DensPart class.")
        elif scheme == "H":
            part = DensPart.from_molecule(molecule, grid=grid, scheme="h", local=False)
        elif scheme == "HI":
            part = DensPart.from_molecule(molecule, grid=grid, scheme="hi", local=False)
        elif scheme is not None and scheme not in ["H", "HI"]:
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

        return cls(molecule, basis, one_rdm, grid, part, grid_2=grid_2, part_2=part_2)

    @classmethod
    def from_file(cls, fname, scheme=None):
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
        grid = get_molecular_grid(molecule, R=1)
        grid_2 = get_molecular_grid(molecule, R=2)

        return cls.from_molecule(molecule, grid, scheme=scheme, grid_2=grid_2)

    def compare_numerical_analytical(self, threshold):
        """
        for the total values of each energy term, calculate analytical and numerical values and check their agreement.
        """

        ### extract analytical terms
        analytical = self.total_analytical()
        nn_total_analytical = analytical["nn_analytical"]
        en_total_analytical = analytical["en_analytical"]
        kin_total_analytical = analytical["kin_analytical"]
        ex_total_analytical = analytical["ex_analytical"]
        coul_total_analytical = analytical["coul_analytical"]

        ### extract numerical terms
        numerical = self.total_numerical()
        nn_total_numerical = numerical["nn_numerical"]
        en_total_numerical = numerical["en_numerical"]
        kin_total_numerical = numerical["kin_numerical"]
        ex_total_numerical = numerical["ex_numerical"]
        coul_total_numerical = numerical["coul_numerical"]
        kin_density = numerical["kin_density"]
        ex_raw = numerical["ex_raw"]
        coul_raw = numerical["coul_raw"]

        ## differences in kcal/mol
        nn_diff = np.abs(nn_total_analytical - nn_total_numerical) / kcalmol
        en_diff = np.abs(en_total_analytical - en_total_numerical) / kcalmol
        kin_diff = np.abs(kin_total_analytical - kin_total_numerical) / kcalmol
        ex_diff = np.abs(ex_total_analytical - ex_total_numerical) / kcalmol
        coul_diff = np.abs(coul_total_analytical - coul_total_numerical) / kcalmol

        ## extract the maximum error, if it's higher than threshold rise error
        max_diff = np.max([nn_diff, en_diff, kin_diff, ex_diff, coul_diff])

        if max_diff > threshold:
            print(
                f"Error between analytical and numertical total nucleus-nucleus repulsion: {nn_diff}"
            )
            print(
                f"Error between analytical and numertical total electron-nucleus attraction: {en_diff}"
            )
            print(f"Error between analytical and numertical total kinetic energy: {kin_diff}")
            print(f"Error between analytical and numertical total exchange energy: {ex_diff}")
            print(f"Error between analytical and numertical total coulumb energe: {coul_diff}")
            raise ValueError(
                f"The difference between analytical and numerical energies exceed the threshold. Maximum allowed: {threshold}, got {max_diff}."
            )

        ### return the integrands to pass for decomposition
        integrands = {
            "kin_density": kin_density,
            "ex_raw": ex_raw,
            "coul_raw": coul_raw,
        }
        totals = {
            "nn_total_analytical": nn_total_analytical,
            "nn_total_numerical": nn_total_numerical,
            "en_total_analytical": en_total_analytical,
            "en_total_numerical": en_total_numerical,
            "kin_total_analytical": kin_total_analytical,
            "kin_total_numerical": kin_total_numerical,
            "ex_total_analytical": ex_total_analytical,
            "ex_total_numerical": ex_total_numerical,
            "coul_total_analytical": coul_total_analytical,
            "coul_total_numerical": coul_total_numerical,
        }
        return integrands, totals

    def total_numerical(self):
        """
        calculate the total numerical value for each term, NN, EN, EE, Kinetic
        """
        ## EN
        rij = np.linalg.norm(self.molecule.coordinates[:, None, :] - self.grid.points, axis=-1)
        en_total_numerical = self.grid.integrate(
            np.sum((-self.molecule.numbers[:, None] * (self.dens / rij)), axis=0)
        )

        ## Kinietic
        kin_density = evaluate_general_kinetic_energy_density(
            self.dm, self.basis, self.grid.points, -0.25
        )
        kin_total_numerical = self.grid.integrate(kin_density)

        ## Coulobm and Exchange
        nao = self.dm.shape[0]
        eval_ao = evaluate_basis(self.basis, self.grid.points)
        # The broadcast operation is the same as the for loop. Summing over mu because broadcasting
        # allocates too much memory
        eval_mo = np.zeros((self.molecule._iodata.mo.coeffs.T.shape[0], eval_ao.shape[1]))
        for i in range(self.molecule._iodata.mo.coeffs.T.shape[0]):
            mo = np.zeros((eval_ao.shape[1]))
            for mu in range(eval_ao.shape[0]):
                mo += self.molecule._iodata.mo.coeffs.T[i, mu] * eval_ao[mu, :]
            eval_mo[i] = mo
        istart = 0
        chunk_size = 250000000 // (nao**2)
        # Establish all ranges for istart/iend
        segments = []
        istart = 0
        i = 0
        while istart < self.grid.points.shape[0]:
            iend = min(istart + chunk_size, self.grid.points.shape[0])
            # segments.append((istart, iend))
            segments.append(
                (i, istart, iend, eval_mo, self.grid, self.dens, self.molecule._iodata, self.dm)
            )
            istart = iend
            i += 1
        # Checking processors
        ncpus = int(os.environ.get("SLURM_CPUS_PER_TASK", default=2))
        # set_start_method("fork")
        ctx = mp.get_context("fork")
        pool = ctx.Pool(processes=ncpus)
        results = []
        # Parallel execution of _integral_point_charge_electron
        for sub in range(0, len(segments), ncpus - 1):
            results_sub_raw = [
                pool.apply_async(_integral_point_charge_electron, args=s)
                for s in segments[sub : sub + (ncpus - 1)]
            ]
            results_sub = [r.get() for r in results_sub_raw]
            results.extend(results_sub)

        pool.close()
        pool.join()

        total_coul_raw = np.zeros(self.grid.points.shape[0])
        total_exch_raw = np.zeros(self.grid.points.shape[0])
        # distribute chunk grid point charges integral outputs
        for out in results:
            total_coul_raw[out[0] : out[1]] = out[2]
            total_exch_raw[out[0] : out[1]] = out[3]

        coul_total_numerical = self.grid.integrate(total_coul_raw)
        ex_total_numerical = self.grid.integrate(total_exch_raw)

        total_numerical = {
            "nn_numerical": np.sum(self.compute_nn_pairwise_array()),
            "en_numerical": en_total_numerical,
            "kin_density": kin_density,
            "kin_numerical": kin_total_numerical,
            "ex_raw": total_exch_raw,
            "ex_numerical": ex_total_numerical,
            "coul_raw": total_coul_raw,
            "coul_numerical": coul_total_numerical,
        }
        return total_numerical

    def total_analytical(self):
        ## NN (this is going to be useless); I agree, we should remove it later
        energy_nn = np.sum(self.compute_nn_pairwise_array())

        ## EN
        en_int = nuclear_electron_attraction_integral(
            self.basis, self.molecule.coordinates, self.molecule.numbers
        )
        en_total_analytical = np.trace(self.dm.dot(en_int))

        ## Kinetic
        kin_int = kinetic_energy_integral(self.basis)
        kin_total_analytical = np.trace(self.dm.dot(kin_int))

        ## Exchange
        ee_int = electron_repulsion_integral(self.basis)
        exch = np.einsum("rs,mnrs->mn", self.dm, ee_int)
        ex_total_analytical = -0.25 * np.trace(self.dm.dot(exch))

        ## Coulomb
        coul = np.einsum("rs,mrns->mn", self.dm, ee_int)
        coul_total_analytical = 0.5 * np.trace(self.dm.dot(coul))

        total_analytical = {
            "nn_analytical": energy_nn,
            "en_analytical": en_total_analytical,
            "kin_analytical": kin_total_analytical,
            "ex_analytical": ex_total_analytical,
            "coul_analytical": coul_total_analytical,
        }
        return total_analytical

    def run_atomic(self, dft_exch=None, dft_corr=None):
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

        # Initialize results dict
        iqa_results = {}

        logging.info("INITIALIZING INTERACTING QUANTUM ATOMS(IQA) ATOMIC CALCULATION")
        iqa_results["nn_total"] = np.sum(self.compute_nn_pairwise_array())
        iqa_results["en_atomic"] = self.compute_en_atomic()
        # (
        #     iqa_results["kin_atomic"],
        #     iqa_results["kin_total_posdef"],
        #     iqa_results["kin_atomic_posdef"],
        # ) = self.kin_iqa()
        iqa_results["kin_atomic"] = self.compute_kin_atomic()
        dft_xc_edens = {}
        # assuming dft_corr and dft_exch specified together
        if dft_corr and dft_exch:
            (
                iqa_results["x_hf_total"],
                iqa_results["coul_total"],
                iqa_results["x_hf_atomic"],
                iqa_results["coul_atomic"],
            ) = self.ee_iqa_hf()
            iqa_results["c_total"], iqa_results["c_atomic"], coeff_mix = self.ee_iqa_dft(
                self.dens, dft_corr
            )
            iqa_results["x_total"], iqa_results["x_atomic"], coeff_mix = self.ee_iqa_dft(
                self.dens, dft_exch
            )
            if "coeff_mix" in dft_xc_edens.keys() and dft_xc_edens["coeff_mix"]:
                if dft_xc_edens["coeff_mix"]:
                    # Scaling HF exchange
                    iqa_results["x_hf_total"] = (
                        iqa_results["x_hf_total"] * dft_xc_edens["coeff_mix"]
                    )
            else:
                iqa_results.pop("x_hf_total")
                iqa_results.pop("x_hf_atomic")

        elif dft_exch:
            (
                iqa_results["x_hf_total"],
                iqa_results["coul_total"],
                iqa_results["x_hf_atomic"],
                iqa_results["coul_atomic"],
            ) = self.ee_iqa_hf()
            iqa_results["xc_total"], iqa_results["xc_atomic"], coeff_mix = self.ee_iqa_dft(
                self.dens, dft_exch
            )
            if coeff_mix:
                # Scaling HF exchange
                iqa_results["x_hf_total"] = iqa_results["x_hf_total"] * coeff_mix
                iqa_results["x_hf_atomic"] = iqa_results["x_hf_atomic"] * coeff_mix
            else:
                iqa_results.pop("x_hf_total")
                iqa_results.pop("x_hf_atomic")
        elif self.molecule.lot == "rhf":
            (
                iqa_results["x_atomic"],
                iqa_results["coul_atomic"],
            ) = self.ee_iqa_hf()
        else:
            raise ValueError(f"Need to specify an exchange functional too. Got {dft_exch}")

        return iqa_results

    def iqa_pairwise(self, dft_exch=None, dft_corr=None, hf_int=False, ee_int=False):
        r"""Compute pairwise interaction for exchange-correlation energy density from DFT
        functionals using BOD partition method.
        """

        logging.info("INITIALIZING INTERACTING QUANTUM ATOMS(IQA) PAIRWISE CALCULATION")
        # assuming dft_corr and dft_exch specified together
        if dft_corr and dft_exch:
            dft_exch_ab, bod = self.ee_iqa_dft_pairwise(dft_exch)
            dft_coul_ab, bod = self.ee_iqa_dft_pairwise(dft_corr)
            logging.info("CALCULATION: DFT EXCHANGE/CORRELATION PAIRWISE")
            print("dft_exch_ab")
            print(dft_exch_ab)
            print("dft_coul_ab")
            print(dft_coul_ab)
            print("dft_total_ab")
            print(dft_exch_ab + dft_coul_ab)
            if hf_int:
                logging.info("CALCULATION: HF COULOMB AND EXCHANGE PAIRWISE")
                logging.warning("6n integrals can be long")
                ab_hf_coul, ab_hf_exch = self.ee_iqa_hf_pairwise()
                print("ab_hf_coul")
                print(ab_hf_coul)
                print("ab_hf_exch")
                print(ab_hf_exch)
                iqa_pairwise_results = {
                    "dft_exch_ab": dft_exch_ab,
                    "dft_coul_ab": dft_coul_ab,
                    "ab_hf_exch": ab_hf_exch,
                    "ab_hf_coul": ab_hf_coul,
                    "bod": bod,
                }
            else:
                logging.info("SKIPPED: HF COULOMB AND EXCHANGE PAIRWISE")
                logging.info("ONLY AB TERMS FOR DFT")
                iqa_pairwise_results = {
                    "dft_exch_ab": dft_exch_ab,
                    "dft_coul_ab": dft_coul_ab,
                    "bod": bod,
                }

            return iqa_pairwise_results

        elif dft_exch:
            dft_xc_ab, bod = self.ee_iqa_dft_pairwise(dft_exch)
            print("dft_xc_ab")
            print(dft_xc_ab)
            if hf_int:
                logging.info("CALCULATION: HF COULOMB AND EXCHANGE PAIRWISE")
                logging.warning("6n integrals can be long")
                ab_hf_coul, ab_hf_exch = self.ee_iqa_hf_pairwise()
                print("ab_hf_coul")
                print(ab_hf_coul)
                print("ab_hf_exch")
                print(ab_hf_exch)
                iqa_pairwise_results = {
                    "dft_xc_ab": dft_xc_ab,
                    "ab_hf_exch": ab_hf_exch,
                    "ab_hf_coul": ab_hf_coul,
                    "bod": bod,
                }
            else:
                logging.info("SKIPPED: HF COULOMB AND EXCHANGE NUM INTEGRATION AB TERMS")
                logging.info("ONLY AB TERMS FOR DFT")
                iqa_pairwise_results = {"dft_xc_ab": dft_xc_ab, "bod": bod}
        else:
            raise ValueError("A libxc must be specified")

        # Print results
        if ee_int:
            from gbasis.integrals.libcint import CBasis

            logging.info("CALCULATION: INTEGRAL EE COULOMB/EXCHANGE")
            flag_libcint = len(
                glob(join(dirname(gbasis.__file__), "integrals", "lib", "libcint.so*"))
            )
            # DEBUG libcint
            # print(flag_libcint)
            if flag_libcint > 0:
                logging.info("USING LIBCINT LIBRARY")
                atoms = [num2sym[z] for z in self.molecule.numbers]
                # check if any cartesian shell
                cart_flag = [s.coord_type for s in self.basis]
                # DEBUG cart/sph
                # print(cart_flag)
                # todo: double check below filter
                if "cartesian" not in cart_flag:
                    cbasis = CBasis(
                        self.basis, atoms, self.molecule.coordinates, coord_type="cartesian"
                    )
                else:
                    cbasis = CBasis(self.basis, atoms, self.molecule.coordinates)

                # Computing ee integrals
                int2e_mo = cbasis.electron_repulsion_integral(
                    transform=self.molecule._iodata.mo.coeffs.T, notation="chemist"
                )
            else:
                logging.info("USING GBASIS PYTHON IMPLEMENTATION")
                logging.warning("EXPECT LONG COMPUTATION")
                int2e_mo = electron_repulsion_integral(
                    self.basis, transform=self.molecule._iodata.mo.coeffs.T, notation="chemist"
                )

            j_coul = 0
            k_ex = 0
            # Mask only occupied Molecular orbitals
            occ_mo = self.molecule._iodata.mo.occs[self.molecule._iodata.mo.occs > 0].shape[0]
            for i in range(occ_mo):
                for j in range(occ_mo):
                    j_coul += 2 * int2e_mo[i, i, j, j]
                    k_ex += int2e_mo[i, j, i, j]

        print("SUMMARY IQA-AB TERMS")

        print("DFT AA/AB TERMS")
        print("DFT ex parition based on:")
        print(
            "Gimferrer M, Salvador P. Exact decompositions of the total KS-DFT exchange–correlation energy into one-and two-center terms. The Journal of Chemical Physics. 2023 Jun 21;158(23)."
        )
        print(
            "Salvador P, Mayer I. One-and two-center physical space partitioning of the energy in the density functional theory. The Journal of chemical physics. 2007 Jun 21;126(23)."
        )
        if dft_exch and dft_corr:
            print("DTF Exchange energy: ", np.sum(iqa_pairwise_results["dft_x_ab"]))
            print("DTF Correlation energy: ", np.sum(iqa_pairwise_results["dft_c_ab"]))
            print()
            print("AA/AB terms DFT Exchange energy:")
            for idx1, at1 in enumerate(self.molecule.numbers):
                for idx2 in range(idx1 + 1):
                    print(
                        f"{at1}({idx1+1}){self.molecule.numbers[idx2]}({idx2+1}) {iqa_pairwise_results['dft_x_ab'][idx1][idx2]}"
                    )
            print("AA/AB terms DFT Correlation energy:")
            for idx1, at1 in enumerate(self.molecule.numbers):
                for idx2 in range(idx1 + 1):
                    print(
                        f"{at1}({idx1+1}){self.molecule.numbers[idx2]}({idx2+1}) {iqa_pairwise_results['dft_c_ab'][idx1][idx2]}"
                    )

        if dft_exch and not dft_corr:
            print("DTF Exchange-Correlation energy: ", np.sum(iqa_pairwise_results["dft_xc_ab"]))
            print()
            print("Atomic DFT Exchange-Correlation energy:")
            print("AA/AB terms DFT Exchange energy:")
            for idx1, at1 in enumerate(self.molecule.numbers):
                for idx2 in range(idx1 + 1):
                    print(
                        f"{at1}({idx1+1}){self.molecule.numbers[idx2]}({idx2+1}) {iqa_pairwise_results['dft_xc_ab'][idx1][idx2]}"
                    )
        print("------------------------------------")
        if hf_int:
            print("Coulomb energy(NUM): ", np.sum(iqa_pairwise_results["ab_hf_coul"]))
            if ee_int:
                print("Coulomb energy(GBASIS): ", j_coul)
                print("DIFF (AU) :", np.sum(iqa_pairwise_results["ab_hf_coul"]) - j_coul)
                print(
                    "DIFF (Kcal) :",
                    (np.sum(iqa_pairwise_results["ab_hf_coul"]) - j_coul) / 0.0015936014376406278,
                )
            print("Exchange energy(NUM): ", np.sum(iqa_pairwise_results["ab_hf_exch"]))
            if ee_int:
                print("HF Exchange energy(GBASIS): ", k_ex)
                print("DIFF (AU) :", np.sum(iqa_pairwise_results["ab_hf_exch"]) - k_ex)
                print(
                    "DIFF (Kcal) :",
                    (np.sum(iqa_pairwise_results["ab_hf_exch"]) - k_ex) / 0.0015936014376406278,
                )

            print()
            print("AA/AB terms Coulomb energy:")
            for idx1, at1 in enumerate(self.molecule.numbers):
                for idx2 in range(idx1 + 1):
                    print(
                        f"{at1}({idx1+1}){self.molecule.numbers[idx2]}({idx2+1}) {iqa_pairwise_results['ab_hf_coul'][idx1][idx2]}"
                    )
            print("AA/AB terms HF exchange energy:")
            for idx1, at1 in enumerate(self.molecule.numbers):
                for idx2 in range(idx1 + 1):
                    print(
                        f"{at1}({idx1+1}){self.molecule.numbers[idx2]}({idx2+1}) {iqa_pairwise_results['ab_hf_exch'][idx1][idx2]}"
                    )
            print()

        return iqa_pairwise_results

    def compute_nn_pairwise_array(self):
        r"""Compute pairwise nuclear-nuclear repulsion energy terms.

        math::
            \sum_{A>B} \frac{Z_{A}Z_{B}}{R_{A}-R_{B}}

        Returns
        -------
        nn_array: np.array(npairs,)
            1D array with nuclear-nuclear repulsion energies for unique pairs of atoms.
            For example, for a molecule with m atoms, the length of the 1D array is m*(m-1)/2.
        """
        logging.info("CALCULATING NUCLEAR-NUCLEAR REPULSION ENERGY TERMS")
        rab = np.triu(
            np.linalg.norm(self.molecule.coordinates[:, None] - self.molecule.coordinates, axis=-1)
        )
        nn_array = np.triu(self.molecule.numbers[:, None] * self.molecule.numbers)[
            np.where(rab > 0)
        ]
        nn_array = nn_array / rab[rab > 0]
        return nn_array

    def compute_kin_atomic(self):
        r"""Decompose kinetic energy into atomic contributions.

        math::
            T_+ (\mathbf{r}_n) = \frac{1}{2} \left. \nabla_{\mathbf{r}}^{2} \gamma(\mathbf{r}, \mathbf{r}')
                                 \right|_{\mathbf{r} = \mathbf{r}' = \mathbf{r}_n}
            Kinetic= T_{\alpha} (\mathbf{r}_n) = T_+(\mathbf{r}_n) + \alpha \nabla^2 \rho(\mathbf{r}_n)

        Where $T_+$ is the positive definite kinetic energy density, with $\gamma$ being the one-electron
        density matrix, and Kinetic is the expression use to return `total_en`.

        Returns
        -------
        kin_atomic: np.array(m, )
            Atomic condensed kinetic energies for `m` atoms in the molecule.
        """
        logging.info("DECOMPOSE KINETIC ENERGY INTO ATOMIC CONTRIBUTIONS")
        # output_posdef = evaluate_posdef_kinetic_energy_density(self.dm, self.basis, self.grid.points)
        # total_kin = self.total_numerical['kin_total']
        # total_kin_posdef = self.grid.integrate(output_posdef)
        if self.part is None:
            raise ValueError("Argument scheme=None, so kinetic energy cannot be composed.")
        kin_atomic = self.part.condense_to_atoms(self.integrands["kin_density"])
        return kin_atomic

    def compute_en_atomic(self, share_factor=0.5):
        r"""Decompose electron-nuclear attraction energy into atomic contributions.

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
        en_atomic: np.array(m, )
            Atomic condensed electron-nuclear attraction energies for `m` atoms in the molecule.
        """
        en_pairwise_matrix = self.compute_en_pairwise_matrix()
        logging.info(f"CALCULATE ELECTRON-NUCLEAR ATOMIC ENERGY TERMS (share={share_factor})")
        # calculate the share matrix for each atom
        share_matrix = np.zeros((self.molecule.natom, self.molecule.natom, self.molecule.natom))
        for i in range(self.molecule.natom):
            share_matrix[i, i, :] = share_factor
            share_matrix[i, :, i] = 1 - share_factor
            share_matrix[i, i, i] = 1
        # compute the share of each pairwise term for each atom
        en_atomic = np.sum(en_pairwise_matrix[None, :, :] * share_matrix, axis=(1, 2))
        return en_atomic

    def compute_en_pairwise_matrix(self):
        """Compute pairwise electron-nuclear attraction energy matrix.

        Returns
        -------
        en_pairwise_matrix: np.array(m, m)
            2D array with pairwise electron-nuclear energies for m atoms in the molecule.
            Diagonal elements correspond to intra-atomic terms (electron density of atom i
            interacting with its own nucleus), while off-diagonal elements correspond to
            inter-atomic terms (electron density of atom i interacting with nucleus of atom j).
        """
        if self.part is None:
            raise ValueError("Argument scheme=None, so electron-nuclear energy cannot be decomposed.")

        logging.info("DECOMPOSE ELECTRON-NUCLEUS ENERGY INTO INTRA- and INTER-ATOMIC TERMS")
        # calculate distance from each atom to each grid point
        rij = np.linalg.norm(self.molecule.coordinates[:, None, :] - self.grid.points, axis=-1)
        # calculate atomic electron-nuclear attraction energy density matrix
        en_pairwise_density = -self.molecule.numbers[None, :, None] * (
            (self.dens[None:,] * self.part.at_weights)[:, None, :] / rij[None, :, :]
        )
        # integrate energy density to get pairwise electron-nuclear energy matrix
        en_pairwise_matrix = np.zeros((self.molecule.natom, self.molecule.natom))
        for i in range(self.molecule.natom):
            for j in range(self.molecule.natom):
                en_pairwise_matrix[i][j] = self.grid.integrate(en_pairwise_density[i][j])

        return en_pairwise_matrix

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

        logging.info("CALCULATING COULOMB AND HF EXCHANGE ENERGY")

        natoms = self.molecule.numbers.shape[0]
        total_coul_raw = self.integrands["coul_raw"]
        total_exch_raw = self.integrands["ex_raw"]

        at_exch = None
        at_coulomb = None
        if self.part:
            if self.part.__class__.__name__ in ["VarHirshfeld", "HirshfeldI", "Hirshfeld"]:
                at_weights = self.part.weights
                at_coulomb_raw = at_weights * total_coul_raw[None, :]
                at_exch_raw = at_weights * total_exch_raw[None, :]
                at_coulomb = np.array(
                    [self.grid.integrate(at_coulomb_raw[i]) for i in range(natoms)]
                )
                at_exch = np.array([self.grid.integrate(at_exch_raw[i]) for i in range(natoms)])
            else:
                at_coulomb = self.part.condense_to_atoms(total_coul_raw)
                at_exch = self.part.condense_to_atoms(total_exch_raw)
            logging.info("Decomposing Coulomb and Exchange into atomic contributions.")
        return at_exch, at_coulomb

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

        # Draft 6N integration
        if not self.grid_2:
            raise ValueError(
                f"Two grids are needed to perform 6N ee integration. Got {self.grid_2}"
            )
        if not self.part_2:
            raise ValueError(
                f"Two partition objects are needed to perform 6N ee integration. Got {self.part_2}"
            )

        points1 = self.grid._grid.points
        points2 = self.grid_2._grid.points
        # Grid1
        eval_ao_g1 = evaluate_basis(self.basis, points1)
        eval_mo_g1 = compute_molecular_orbitals_from_ao(self.molecule, eval_ao_g1)
        # Grid2
        eval_ao_g2 = evaluate_basis(self.basis, points2)
        eval_mo_g2 = compute_molecular_orbitals_from_ao(self.molecule, eval_ao_g2)

        # Compute density/exchange density from Molecular orbitals
        # Because only restricted _occs_a.shape[0] == _occs_b.shape[0]
        occupied_mo = self.molecule._iodata.mo.occs[self.molecule._iodata.mo.occs > 0].shape[0]

        if self.part:
            if self.part_2:
                logging.warning("Molecular 6N integration using local grids")
                ab_hf_coul = np.zeros((len(self.molecule.numbers), len(self.molecule.numbers)))
                ab_hf_exch = np.zeros((len(self.molecule.numbers), len(self.molecule.numbers)))
                natoms = len(self.molecule.numbers)
            else:
                raise ValueError("6N integration uses local grids and 2 partition objects")

        for a in range(natoms):
            print(a)
            atgrid_a = self.grid._grid.get_localgrid(self.molecule.coordinates[a], 300)
            for b in range(a + 1):
                print(b)
                atgrid_b = self.grid_2.get_localgrid(self.molecule.coordinates[b], 300)

                integral_coul_ab = 0
                integral_ex_ab = 0
                for id1, p1 in enumerate(atgrid_a.indices):
                    progress = p1 / atgrid_a.points.shape[0]
                    if progress * 100 in [15.0, 25.0, 50.0, 75.0, 90.0, 100.0]:
                        logging.info(f"PROGRESS 6N INTEGRATION: {progress * 100}%")
                    d_coul = 0
                    d_ex = 0
                    for i in range(occupied_mo):
                        for j in range(occupied_mo):
                            d_coul += (
                                2
                                * eval_mo_g1[i, p1]
                                * eval_mo_g1[i, p1]
                                * eval_mo_g2[j]
                                * eval_mo_g2[j]
                            )
                            d_ex += (
                                eval_mo_g1[i, p1]
                                * eval_mo_g1[j, p1]
                                * eval_mo_g2[j]
                                * eval_mo_g2[i]
                            )

                    rij = np.linalg.norm(self.grid.points[p1, :] - points2, axis=-1)
                    rij[rij == 0] = 1.0e-9
                    d_ex = d_ex / rij
                    d_coul = d_coul / rij
                    # Subset part weights b
                    w_subset_b = self.part_2.at_weights[b][atgrid_b.indices]
                    # Subset d_coul and d_ex
                    d_coul_b = d_coul[atgrid_b.indices]
                    d_ex_b = d_ex[atgrid_b.indices]
                    # Atomic integration for r2 using part weights
                    part_coul_b = atgrid_b.integrate(d_coul_b * w_subset_b)
                    part_ex_b = atgrid_b.integrate(d_ex_b * w_subset_b)
                    integral_coul_ab += part_coul_b * (
                        atgrid_a.weights[id1] * self.part.at_weights[a][p1]
                    )
                    integral_ex_ab += part_ex_b * (
                        atgrid_a.weights[id1] * self.part.at_weights[a][p1]
                    )

                ab_hf_coul[a, b] = integral_coul_ab
                ab_hf_exch[a, b] = integral_ex_ab
                if a != b:
                    ab_hf_coul[b, a] = integral_coul_ab
                    ab_hf_exch[b, a] = integral_ex_ab

        return ab_hf_coul, ab_hf_exch

    def ee_iqa_dft(self, dft_dens, libxc_label):

        natoms = self.molecule.numbers.shape[0]

        # Test pylibxc
        print(libxc_label)
        func = pylibxc.LibXCFunctional(libxc_label, "unpolarized")

        input_libxc = {}
        # From pylibxc flags.py
        # XC_FAMILY_LDA = 1
        # XC_FAMILY_HYB_LDA = 128
        if func.get_family() in [1, 128]:
            input_libxc["rho"] = self.dens
        # From pylibxc flags.py
        # XC_FAMILY_GGA = 2
        # XC_FAMILY_HYB_GGA = 32
        elif func.get_family() in [2, 32]:
            # Getting the norm-squared of the gradient of the spin-summed electron density(sigma_full)
            grad = evaluate_density_gradient(self.dm, self.basis, self.grid.points)
            sigma_full = 4 * ((grad / 2) ** 2).sum(axis=1)
            input_libxc["rho"] = self.dens
            input_libxc["sigma"] = sigma_full

        results_func = func.compute(input_libxc)
        dft_xc_dens = results_func["zk"].reshape(-1)

        # From pylibxc flags.py
        # XC_FAMILY_HYB_GGA = 32
        # XC_FAMILY_HYB_MGGA = 64
        # XC_FAMILY_HYB_LDA = 128
        if func.get_family() in [32, 128]:
            coeff_mix = func.get_hyb_exx_coef()
        else:
            coeff_mix = None

        dft_xc_total = self.grid.integrate(self.dens, dft_xc_dens)

        at_dft_xc = None
        if self.part:
            if self.part.__class__.__name__ in ["VarHirshfeld", "HirshfeldI", "Hirshfeld"]:
                at_weights = self.part.weights
                at_dft_xc_raw = at_weights * dft_xc_dens[None, :]
                at_dft_xc = np.array(
                    [self.grid.integrate(at_dft_xc_raw[i], self.dens) for i in range(natoms)]
                )
            else:
                at_dft_xc = self.part.condense_to_atoms((dft_xc_dens * self.dens))
            logging.info("Decomposing XC into atomic contributions.")
            print("XC component")
            print(at_dft_xc)
            print(np.sum(at_dft_xc), dft_xc_total)
            assert_almost_equal(np.sum(at_dft_xc), dft_xc_total, decimal=2)

        print()
        return dft_xc_total, at_dft_xc, coeff_mix

    def ee_iqa_dft_pairwise(self, libxc_label):

        # Get different data
        at_weights = self.part.weights
        natoms = self.molecule.numbers.shape[0]

        # Get functional form pylibxc
        func = pylibxc.LibXCFunctional(libxc_label, "unpolarized")

        # BOD for LDA
        # Evaluate Atomic orbitals basis
        eval_ao = evaluate_basis(self.basis, self.grid.points)
        # Get Molecular orbitals from AO
        eval_mo = compute_molecular_orbitals_from_ao(self.molecule, eval_ao)
        # Computing Molecular overlap from eval_mo
        sij = eval_mo[None, :, :] * eval_mo[:, None, :]
        # Getting atomic overlap matrices
        sij_at_1 = at_weights[:, None, None, :] * sij[None, :, :, :]

        # Computing BOD for each pair of atoms
        bod = np.zeros((natoms, natoms, self.grid.points.shape[0]))
        for a in range(natoms):
            for b in range(natoms):
                for i in range(
                    self.molecule._iodata.mo.occs[self.molecule._iodata.mo.occs > 0].shape[0]
                ):
                    for j in range(
                        self.molecule._iodata.mo.occs[self.molecule._iodata.mo.occs > 0].shape[0]
                    ):
                        bod[a, b] += (
                            2
                            * (
                                at_weights[a] * self.grid.integrate(sij_at_1[b][i][j])
                                + at_weights[b] * self.grid.integrate(sij_at_1[a][j][i])
                            )
                            * eval_mo[i]
                            * eval_mo[j]
                        )
        # Compute Grad BOD for GGA and Hyb
        # From pylibxc flags.py
        # XC_FAMILY_GGA = 2
        # XC_FAMILY_HYB_GGA = 32
        if func.get_family() in [2, 32]:
            # BOD for GGA and Hybrid GGA
            # Evaluate molecular orbitals gradient
            # todo: this could potentially go to utils
            orders_d = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            eval_mo_grad = []
            eval_mo_grad_condensed = []
            for ind, orders in enumerate(orders_d):
                eval_mo_grad_order = evaluate_deriv_basis(
                    self.basis,
                    self.grid.points,
                    np.array(orders),
                    transform=self.molecule._iodata.mo.coeffs.T,
                )
                eval_mo_grad.append(eval_mo_grad_order[:, :, None])
                eval_mo_grad_condensed.append(eval_mo_grad_order)
            eval_mo_grad = np.concatenate(eval_mo_grad, axis=2)

            # Computing BOD for GGA functionals for each pair of atoms
            print(natoms, natoms, self.grid.points.shape[0])
            grad_bod = np.zeros((natoms, natoms, self.grid.points.shape[0]))
            for a in range(natoms):
                for b in range(natoms):
                    for i in range(
                        self.molecule._iodata.mo.occs[self.molecule._iodata.mo.occs > 0].shape[0]
                    ):
                        for j in range(
                            self.molecule._iodata.mo.occs[self.molecule._iodata.mo.occs > 0].shape[
                                0
                            ]
                        ):
                            grad_bod[a, b] += np.sum(
                                (
                                    at_weights[a, :, None] * self.grid.integrate(sij_at_1[b][i][j])
                                    + at_weights[b, :, None]
                                    * self.grid.integrate(sij_at_1[a][i][j])
                                )
                                * eval_mo_grad[i, :, :]
                                * eval_mo[j, :, None],
                                axis=-1,
                            )

            grad_bod = 16 * (grad_bod**2)

        # Integration BOD
        int_bod = 0
        for a in range(natoms):
            for b in range(a):
                int_bod += self.grid.integrate(bod[a, b])
        print("Integration BOD (TOTAL)", int_bod)

        # Combine atomic partition weights
        combine_atweights = np.zeros((at_weights[0].shape[0]))
        # Below code for grid is for Horton grid.
        for at in range(natoms):
            start = self.grid._grid.indices[at]
            end = self.grid._grid.indices[at + 1]
            combine_atweights[start:end] = at_weights[at][start:end]

        # LDA
        input_libxc = {}
        # From pylibxc flags.py
        # XC_FAMILY_LDA = 1
        # XC_FAMILY_HYB_LDA = 128
        if func.get_family() in [1, 128]:
            dft_pairwise = np.zeros((natoms, natoms))
            for a in range(natoms):
                for b in range(natoms):
                    input_libxc["rho"] = bod[a, b]
                    results_func = func.compute(input_libxc)
                    dft_xc_dens = results_func["zk"].reshape(-1)
                    dft_pairwise[a, b] = self.grid.integrate(
                        (dft_xc_dens * bod[a, b] * combine_atweights)
                    )

        # GGA and Hyb
        # From pylibxc flags.py
        # XC_FAMILY_GGA = 2
        # XC_FAMILY_HYB_GGA = 32
        elif func.get_family() in [2, 32]:
            dft_pairwise = np.zeros((natoms, natoms))
            for a in range(natoms):
                for b in range(natoms):
                    input_libxc["rho"] = bod[a, b]
                    input_libxc["sigma"] = grad_bod[a, b]
                    results_func = func.compute(input_libxc)
                    dft_xc_dens = results_func["zk"].reshape(-1)
                    dft_pairwise[a, b] = self.grid.integrate(
                        (dft_xc_dens * bod[a, b] * combine_atweights)
                    )

        return dft_pairwise, bod
