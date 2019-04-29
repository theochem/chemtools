"""Class for obtaining quasiatomic basis set transformation and convenient attribute/properties.

This code has been copied from https://github.com/QuantumElephant/dumbo.

"""
import numpy as np
from . import orthogonalization as orth


class QuasiTransformationError(Exception):
    """Error in obtaining the quasi basis set."""

    pass


class QuasiTransformation(object):
    """Container for the methods that generate quasi/pseudo atomic orbitals.

    Since these methods use the same parameters, they were grouped into a class
    as they share attributes.

    Attributes
    ----------
    _coeff_ab_mo : np.ndarray
        Transformation matrix from atomic basis to molecular orbitals
    _olp_ab_ab : np.ndarray(N,N)
        Overlap of the atomic basis functions
    _olp_aao_ab : np.ndarray(N,M)
        Overlap of the the accurate atomic orbitals with the atomic basis functions
    _olp_aao_aao : np.ndarray(M,M)
        Overlap of the accurate atomic orbitals

    Properties
    ----------
    coeff_ab_mo
        Coefficient matrix from atomic basis to molecular orbitals
    olp_ab_ab
        Overlap matrix of the atomic basis functions
    olp_aao_ab
        Overlap matrix of the AAO's with the atomic basis functions
    olp_aao_aao
        Overlap matrix of the AAO functions
    indices_span
        Boolean numpy indices of the MO's that are exactly described by the
        quasi orbitals
    olp_ab_aao
        Matrix of the overlaps between atomic basis functions and accurate
        atomic basis functions
    olp_mo_mo
        Matrix of the overlaps between molecular orbitals
    olp_aao_mo
        Matrix of the overlaps between accurate atomic orbitals and molecular
        orbitals
    olp_mo_aao
        Matrix of the overlaps between molecular orbitals and accurate atomic
        orbitals
    coeff_ab_omo
        Transformation matrix from the atomic basis to the occupied molecular
        orbitals
    olp_omo_omo
        Matrix of the overlaps between the occupied molecular orbitals
    olp_omo_ab
        Matrix of the overlaps between the occupied molecular orbitals and
        the atomic basis functions
    olp_ab_omo
        Matrix of the overlaps between the atomic basis functions and the
        occupied molecular orbitals
    olp_omo_aao
        Matrix of the overlaps between the occupied molecular orbitals and
        the accurate atomic orbitals
    olp_aao_omo
        Matrix of the overlaps between the accurate atomic orbitals and the
        occupied molecular orbitals
    num_ab
        Number of atomic basis functions
    num_aao
        Number of AAO functions
    num_mo
        Number of molecular orbitals

    Methods
    -------
    make_mmo
        Makes a basis set (mMO) that spans the same space as quambo or quao
    quambo
        Gets the transformation from AB to QUAMBO
    quao
        Gets the transformation from AB to QUAO
    iao
        Gets the transformation from AB to IAO

    Notes
    -----
    1. The following notations are used:
        ab = atomic basis
        aao = accurate atomic orbitals (from accurate calculations of an atom)
        mo = molecular orbital
    2. No assumptions are made for the spin of the molecular orbitals. They are
       not assumed to be spatial or spin orbitals. If they are spin orbitals and
       you'd like spin separated quasi orbitals, then the quasi code can be run
       twice, one for the alphas and one for the betas, or a block diagonal
       overlaps should be created (there are alpha and beta atomic basis and
       AAO functions). If they are spin orbitals and you'd like spatial quasi
       orbitals, then you should vertically stack the alpha and the beta blocks
       (NEEDS WORK!!)
    """

    def __init__(self, coeff_ab_mo, olp_ab_ab, olp_aao_ab, olp_aao_aao, indices_span):
        """Initialize attributes.

        Parameters
        ----------
        coeff_ab_mo : np.ndarray(K,N)
            Transformation matrix from atomic basis to molecular orbitals
        olp_ab_ab : np.ndarray(K,K)
            Matrix of the overlaps between atomic basis functions
        olp_aao_ab : np.ndarray(L,K)
            Matrix of the overlaps between accurate atomic orbitals and atomic
            basis functions
        olp_aao_aao : np.ndarray(L,L)
            Matrix of the overlaps between accurate atomic orbitals
        indices_span : np.ndarray(M,)
            Indices that indicate which of the mo's will be exactly described

        Raises
        ------
        TypeError
            If coeff_ab_mo is two dimensional numpy array
            If olp_ab_ab is not a two dimensional square numpy array
            If olp_aao_ab is not a two dimensional numpy array
            If olp_aao_aao is not a two dimensional square numpy array
        ValueError
            If olp_ab_ab is not normalized
            If olp_ab_ab is not symmetric
            If olp_aao_aao is not normalized
            If olp_aao_aao is not symmetric
            If number of ab is not consistent for coeff_ab_mo, olp_ab_ab and olp_aao_ab
            If number of aao is not consistent for olp_aao_ab and olp_aao_aao
            If molecular orbitals are not normalized
        """
        # indices_span is not included because it has a setter
        if not (isinstance(coeff_ab_mo, np.ndarray) and len(coeff_ab_mo.shape) == 2):
            raise TypeError("Given coefficient matrix is not a two dimensional numpy array")
        if not (
            isinstance(olp_ab_ab, np.ndarray)
            and len(olp_ab_ab.shape) == 2
            and olp_ab_ab.shape[0] == olp_ab_ab.shape[1]
        ):
            raise TypeError(
                "Given overlap matrix for atomic basis is not a two dimensional square "
                "numpy array"
            )
        if not np.allclose(np.diag(olp_ab_ab), np.ones(olp_ab_ab.shape[0])):
            raise ValueError("Given overlap matrix for atomic basis is not normalized")
        if not np.allclose(olp_ab_ab, olp_ab_ab.T):
            raise ValueError("Given overlap matrix for atomic basis is not symmetric")
        if not np.all(orth.eigh(olp_ab_ab)[0] >= 0):
            raise ValueError("Given overlap matrix for atomic basis is not positive semidefinite")

        if not (isinstance(olp_aao_ab, np.ndarray) and len(olp_aao_ab.shape) == 2):
            raise TypeError(
                "Given overlap matrix for atomic basis and AAO is not a two "
                "dimensional numpy array"
            )
        if not (
            isinstance(olp_aao_aao, np.ndarray)
            and len(olp_aao_aao.shape) == 2
            and olp_aao_aao.shape[0] == olp_aao_aao.shape[1]
        ):
            raise TypeError(
                "Given overlap matrix for AAO is not a two dimensional square numpy " "array"
            )
        if not np.allclose(np.diag(olp_aao_aao), np.ones(olp_aao_aao.shape[0])):
            raise ValueError("Given overlap matrix for AAO is not normalized")
        if not np.allclose(olp_aao_aao, olp_aao_aao.T):
            raise ValueError("Given overlap matrix for AAO is not symmetric")
        if not np.all(orth.eigh(olp_aao_aao)[0] >= 0):
            raise ValueError("Given overlap matrix for AAO is not positive semidefinite")
        if not (coeff_ab_mo.shape[0] == olp_ab_ab.shape[0] == olp_aao_ab.shape[1]):
            raise ValueError("Number of atomic basis functions is not consistent")
        if olp_aao_ab.shape[0] != olp_aao_aao.shape[0]:
            raise ValueError("Number of AAO basis functions is not consistent")

        olp_mo_mo = coeff_ab_mo.T.dot(olp_ab_ab).dot(coeff_ab_mo)
        if not np.allclose(np.diag(olp_mo_mo), np.ones(olp_mo_mo.shape[0])):
            raise ValueError(
                "The overlap of the molecular orbitals, calculated from coeff_ab_mo "
                "and olp_ab_ab is not normalized"
            )

        # Initialize private variables
        self._coeff_ab_mo = coeff_ab_mo
        self._olp_ab_ab = olp_ab_ab
        self._olp_aao_ab = olp_aao_ab
        self._olp_aao_aao = olp_aao_aao
        self._concatenate_ab_aao = False
        # Use setter
        self._indices_span = None
        self.indices_span = indices_span

    @property
    def coeff_ab_mo(self):
        """Return the oefficient matrix from atomic basis to molecular orbitals."""
        if not self._concatenate_ab_aao:
            return self._coeff_ab_mo
        else:
            return np.vstack((self._coeff_ab_mo, np.zeros((self.num_aao, self.num_mo))))

    @property
    def olp_ab_ab(self):
        """Return the overlap matrix of the atomic basis functions."""
        if not self._concatenate_ab_aao:
            return self._olp_ab_ab
        else:
            olp_ab_new = np.hstack((self._olp_ab_ab, self._olp_aao_ab.T))
            olp_aao_new = np.hstack((self._olp_aao_ab, self._olp_aao_aao))
            return np.vstack((olp_ab_new, olp_aao_new))

    @property
    def olp_aao_ab(self):
        """Return the overlap matrix of the AAO's with the atomic basis functions."""
        if not self._concatenate_ab_aao:
            return self._olp_aao_ab
        else:
            return np.hstack((self._olp_aao_ab, self._olp_aao_aao))

    @property
    def olp_aao_aao(self):
        """Return the overlap matrix of the AAO functions."""
        return self._olp_aao_aao

    @property
    def indices_span(self):
        """Return the boolean indices of the MO's that are exactly described by the quasi orbitals.
        """
        return self._indices_span

    @indices_span.setter
    def indices_span(self, indices):
        """Set the indices_span.

        Parameters
        ----------
        Indices : np.ndarray of {bool, int}
            The indices of the molecular orbitals that will be exactly described
            by the quasi orbitals

        Raises
        ------
        TypeError
            If indices is not a one dimensional numpy array
            If indices are not boolean or integers
        ValueError
            If indices are integral and outside of valid range (between 0 and num_mo)
            If indices are integral and are repeated
            If indices are boolean and does not have exactly the right size (num_mo)
        """
        if not (isinstance(indices, np.ndarray) and len(indices.shape) == 1):
            raise TypeError("Indices are not given as a one dimensional numpy array")
        if indices.dtype not in [int, bool]:
            raise TypeError("Indices can only be given as a integral or boolean index")
        if indices.dtype == int:
            if not (np.all(indices < self.num_mo) and np.all(indices >= 0)):
                raise ValueError("Invalid index (greater than number of MOs or less than 0))")
            if not len(set(indices)) == indices.size:
                raise ValueError("Given indices contains repeated indices")
            temp = np.zeros(self.num_mo, dtype=bool)
            temp[indices] = True
            indices = temp
        elif indices.dtype == bool:
            if indices.size != self.num_mo:
                raise ValueError(
                    "Given boolean numpy indices needs to have the same number of"
                    " entries as the number of molecular orbitals"
                )
        self._indices_span = indices

    @property
    def olp_ab_aao(self):
        """Return the matrix of the overlap between atomic basis functions (rows) and AAO's (cols).
        """
        return self.olp_aao_ab.T

    @property
    def olp_mo_mo(self):
        """Return the matrix of the overlaps between molecular orbitals."""
        return self.coeff_ab_mo.T.dot(self.olp_ab_ab).dot(self.coeff_ab_mo)

    @property
    def olp_aao_mo(self):
        """Return the matrix of the overlaps between AAO's and molecular orbitals."""
        return self.olp_aao_ab.dot(self.coeff_ab_mo)

    @property
    def olp_mo_aao(self):
        """Return the matrix of the overlaps between molecular orbitals and AAO's."""
        return self.olp_aao_mo.T

    @property
    def coeff_ab_omo(self):
        """Return the transformation matrix from the atomic basis to the occupied MO's."""
        return self.coeff_ab_mo[:, self.indices_span]

    @property
    def olp_omo_omo(self):
        """Return the matrix of the overlaps between the occupied molecular orbitals."""
        return self.coeff_ab_omo.T.dot(self.olp_ab_ab).dot(self.coeff_ab_omo)

    @property
    def olp_omo_ab(self):
        """Return the matrix of the overlaps between the occupied MO's and the AB's."""
        return self.coeff_ab_omo.T.dot(self.olp_ab_ab)

    @property
    def olp_ab_omo(self):
        """Return the matrix of the overlaps between the AB's and the occupied MO's."""
        return self.olp_omo_ab.T

    @property
    def olp_omo_aao(self):
        """Return the matrix of the overlaps between the occupied MO's and AAO's."""
        return self.olp_mo_aao[self.indices_span]

    @property
    def olp_aao_omo(self):
        """Return the matrix of the overlaps between the AAO's and the occupied MO's."""
        return self.olp_omo_aao.T

    @property
    def num_ab(self):
        """Return the number of atomic basis functions."""
        assert self.olp_ab_ab.shape[0] == self.olp_aao_ab.shape[1]
        return self.olp_aao_ab.shape[1]

    @property
    def num_aao(self):
        """Return the number of AAO functions."""
        return self._olp_aao_ab.shape[0]

    @property
    def num_mo(self):
        """Return the number of molecular orbitals."""
        return self._coeff_ab_mo.shape[1]

    @property
    def num_omo(self):
        """Return the number of selected molecular orbitals."""
        return sum(self.indices_span)

    @property
    def num_vmo(self):
        """Return the number of molecular orbitals that were not selected."""
        return sum(-self.indices_span)

    def make_mmo(self, dim_mmo=None, olp_aao_mo=None):
        """Return minimal molecular orbitals.

        Parameters
        ----------
        dim_mmo : {None, int}, optional
            Total dimension of the MMO space (Default is the dimension of the AAO space)
        olp_aao_mo : {None, (N,M) np.ndarray, BlockDiagMatrix}, optional
            Overlap of AAO's with MO's
            Used when the mMO construction uses overlap that is different from
            the one stored as attribute, as is the case for QUAO construction

        Returns
        -------
        coeff_ab_mmo : np.ndarray
            Transformation matrix from atomic basis functions to mMO's

        Raises
        ------
        QuasiTransformationError
            If the dimension of the MMO space is larger than the number of molecular orbitals
            If the dimension of the MMO space is smaller than the space that needs to be spanned
            If the `olp_aao_mo` is not a numpy array

        References
        ----------
        .. [1] Lu. W.C.; Wang, C.Z.; Schmidt, W.; Bytautas, L.;Ho K.M.;
           Ruedenberg, K. Ruedenberg. Molecule intrinsic minimal basis sets. I.
           Exact resolution of ab initio optimized molecular orbitals in terms
           of deformed atomic minimal
        .. [2] West, A.C.; Schmidt, M.W.
           comprehensive analysis of molecule-intrinsic quasiatomic, bonding,
           and correlating orbitals. I. Hartree-Fock wave functions. J. Chem.
           Phys. 2013, 139, 234107.
        .. [3] David's paper
        """
        # Check dim_mmo
        if dim_mmo is None:
            dim_mmo = self.num_aao
        elif dim_mmo > self.num_mo:
            raise QuasiTransformationError(
                "Dimension of MMO space, {0}, is larger than the number "
                "of molecular orbitals, {1} ".format(dim_mmo, self.num_mo)
            )
        elif dim_mmo < np.sum(self.indices_span):
            raise QuasiTransformationError(
                "Dimension of MMO space, {0}, is smaller than the space "
                "you want to span, "
                "{1}".format(dim_mmo, np.sum(self.indices_span))
            )

        # Check olp_aao_mo
        if olp_aao_mo is None:
            olp_aao_mo = self.olp_aao_mo
        elif not isinstance(olp_aao_mo, np.ndarray):
            raise QuasiTransformationError("Given olp_aao_mo is not a numpy array")

        # Set local variables
        dim_span = np.sum(self.indices_span)
        num_to_add = dim_mmo - dim_span

        # Create occupied MMO
        coeff_occmo_occmmo = np.identity(dim_span)

        # Create virtual MMO
        #  find overlap between aao and virtuals
        olp_aao_virmo = olp_aao_mo[:, -self.indices_span]
        #  from the right singular vector of olp_aao_virmo
        temp_right_sing_vec = orth.svd(olp_aao_virmo)[2]
        #  select vectors with largest (num_to_add) singular values
        coeff_virmo_virmmo = temp_right_sing_vec[:num_to_add, :].T

        # Combine the occupied and virtual MMO's
        coeff_mo_occmmo = np.zeros((self.indices_span.size, coeff_occmo_occmmo.shape[1]))
        coeff_mo_occmmo[self.indices_span] = coeff_occmo_occmmo
        coeff_mo_virmmo = np.zeros((self.indices_span.size, coeff_virmo_virmmo.shape[1]))
        coeff_mo_virmmo[-self.indices_span] = coeff_virmo_virmmo
        coeff_mo_mmo = np.hstack((coeff_mo_occmmo, coeff_mo_virmmo))

        # Express MMO wrt atomic basis functions
        return self.coeff_ab_mo.dot(coeff_mo_mmo)

    def quambo(self, dim=None, is_orth=False, basis="ab"):
        """Construct QUAMBO basis.

        Parameters
        ----------
        dim : {None, int}, optional
            Number of QUAMBO basis functions (Default is the number of AAO functions)
        is_orth : bool, optional
            Flag for orthogonalizing quambo orbitals
        basis : str
            Basis set with respect to which the quambo will be described
            'ab' will result in quambo with respect to atomic basis functions
            'mo' will result in quambo with respect to molecular orbitals

        Returns
        -------
        coeff_ab_quambo : np.ndarray
            Transformation matrix from atomic basis functions to QUAMBO's

        Raises
        ------
        QuasiTransformationError
            If the QUAMBO's cannot be orthogonalized.

        Notes
        -----
        .. math::

            \big| \mathrm{QUAMBO}_k \big> &= \mathrm{proj}_{\mathrm{mMO}^\mathrm{QUAMBO}}\\
            &= \sum_i \big| \mathrm{mMO}^\mathrm{QUAMBO}_i \big>
               (olp_{{\mathrm{mMO}^\mathrm{QUAMBO}}}^{-1})_{ij}
               \big< \mathrm{mMO}^\mathrm{QUAMBO}_j \big| \mathrm{AAO}_k \big>\\
            &= \sum_j \big| AB_j \big> (C^{\mathrm{AB},\mathrm{mMO}^\mathrm{QUAMBO}}
               olp_{{\mathrm{mMO}^\mathrm{QUAMBO}}}^{-1}
               olp_{\mathrm{mMO}^\mathrm{QUAMBO}, \mathrm{AAO}})_{ji}\\

        References
        ----------
        .. [1] Lu. W.C.; Wang, C.Z.; Schmidt, W.; Bytautas, L.;Ho K.M.;
           Ruedenberg, K. Ruedenberg. Molecule intrinsic minimal basis sets. I.
           Exact resolution of ab initio optimized molecular orbitals in terms
           of deformed atomic minimal
        """
        # Find MMO for QUAMBOs
        coeff_ab_mmo = self.make_mmo(dim_mmo=dim)

        # Get transformation
        olp_mmo_mmo = coeff_ab_mmo.T.dot(self.olp_ab_ab).dot(coeff_ab_mmo)
        olp_mmo_aao = coeff_ab_mmo.T.dot(self.olp_ab_aao)
        coeff_mmo_proj = project(olp_mmo_mmo, olp_mmo_aao)

        # Orthogonalize
        if is_orth:
            olp_proj_proj = coeff_mmo_proj.T.dot(olp_mmo_mmo).dot(coeff_mmo_proj)
            coeff_mmo_proj = coeff_mmo_proj.dot(orth.power_symmetric(olp_proj_proj, -0.5))
            if not np.allclose(
                np.diag(coeff_mmo_proj.T.dot(olp_mmo_mmo).dot(coeff_mmo_proj)), 1, atol=1e-2
            ):
                raise QuasiTransformationError("QUAMBO's could not be orthogonalized.")
        # Normalize
        olp_proj_proj = coeff_mmo_proj.T.dot(olp_mmo_mmo).dot(coeff_mmo_proj)
        coeff_mmo_proj *= np.diag(olp_proj_proj) ** (-0.5)

        if basis == "ab":
            return coeff_ab_mmo.dot(coeff_mmo_proj)
        elif basis == "mo":
            coeff_mo_ab = project(self.olp_mo_mo, self.coeff_ab_mo.T.dot(self.olp_ab_ab))
            return coeff_mo_ab.dot(coeff_ab_mmo)
        else:
            raise ValueError("`basis` must be one of 'ab' or 'mo'")

    def quao(self, dim=None, is_orth=False, basis="ab"):
        """Construct QUAO basis.

        Parameters
        ----------
        dim : {None, int}, optional
            Number of QUAO basis functions (Default is the number of AAO functions)
        is_orth : bool, optional
            Flag for orthogonalizing quao orbitals
        basis : str
            Basis set with respect to which the quambo will be described
            'ab' will result in quambo with respect to atomic basis functions
            'mo' will result in quambo with respect to molecular orbitals

        Returns
        -------
        coeff_ab_quao
            Transformation matrix from atomic basis functions to QUAO's

        Raises
        ------
        QuasiTransformationError
            If the QUAMBO's cannot be orthogonalized.

        Notes
        -----
        .. math::

            \big| \mathrm{QUAO}_k \big> &= \mathrm{proj}_{\mathrm{mMO}^\mathrm{QUAO}}\\
            &= \sum_i \big| \mathrm{mMO}^\mathrm{QUAO}_i \big>
               (olp_{{\mathrm{mMO}^\mathrm{QUAO}}}^{-1})_{ij}
               \big< \mathrm{mMO}^\mathrm{QUAO}_j \big| \mathrm{AAO}_k \big>\\
            &= \sum_j \big| AB_j \big> (C^{\mathrm{AB},\mathrm{mMO}^\mathrm{QUAO}}
               olp_{{\mathrm{mMO}^\mathrm{QUAO}}}^{-1}
               olp_{\mathrm{mMO}^\mathrm{QUAO}, \mathrm{AAO}})_{ji}\\

        References
        ----------
        .. [1] West, A.C.; Schmidt, M.W. Gordon, M.S; Ruedenberg, K. A
           comprehensive analysis of molecule-intrinsic quasiatomic, bonding,
           and correlating orbitals. I. Hartree-Fock wave functions. J. Chem.
           Phys. 2013, 139, 234107.
        """
        # Orthogonalize AAOs
        olp_oaao_mo = orth.power_symmetric(self.olp_aao_aao, -0.5).dot(self.olp_aao_mo)

        # Get MMOs using the orthogonalized AAOs (MMO for QUAOs)
        coeff_ab_mmo = self.make_mmo(dim_mmo=dim, olp_aao_mo=olp_oaao_mo)

        # Find transformation for QUAOs
        olp_mmo_mmo = coeff_ab_mmo.T.dot(self.olp_ab_ab).dot(coeff_ab_mmo)
        olp_mmo_aao = coeff_ab_mmo.T.dot(self.olp_ab_aao)
        coeff_mmo_proj = project(olp_mmo_mmo, olp_mmo_aao)

        # Orthogonalize
        if is_orth:
            olp_proj_proj = coeff_mmo_proj.T.dot(olp_mmo_mmo).dot(coeff_mmo_proj)
            coeff_mmo_proj = coeff_mmo_proj.dot(orth.power_symmetric(olp_proj_proj, -0.5))
            if not np.allclose(
                np.diag(coeff_mmo_proj.T.dot(olp_mmo_mmo).dot(coeff_mmo_proj)), 1, atol=1e-2
            ):
                raise QuasiTransformationError("QUAMBO's could not be orthogonalized.")
        # Normalize
        olp_proj_proj = coeff_mmo_proj.T.dot(olp_mmo_mmo).dot(coeff_mmo_proj)
        coeff_mmo_proj *= np.diag(olp_proj_proj) ** (-0.5)
        if basis == "ab":
            return coeff_ab_mmo.dot(coeff_mmo_proj)
        elif basis == "mo":
            coeff_mo_ab = project(self.olp_mo_mo, self.coeff_ab_mo.T.dot(self.olp_ab_ab))
            return coeff_mo_ab.dot(coeff_ab_mmo)
        else:
            raise ValueError("`basis` must be one of 'ab' or 'mo'")

    def iao(self, iaotype=1, is_orth=True):
        """Construct the IAO basis.

        Parameters
        ----------
        iaotype : {1, 2}, optional
            If 1 then IAO of form (type 1) is returned
            If 2 then IAO of form (type 2) is returned
            See Notes for the equation

        Returns
        -------
        coeff_basis_iao : np.ndarray
            Transformation matrix from a basis set to IAO's
            If `iaotype` is 1, then basis set is the AB.
            If `iaotype` is 2, then basis set is the concatenation of AB's and AAO's.

        Raises
        ------
        QuasiTransformationError
            If the `iaotype` is not 1 or 2
            If IAO's cannot be orthonormalized

        Notes
        -----
        projections:
        .. math::

            \proj_\ab \ket{\aao_i}
                &= \sum_{kl} \ket{\ab_k} (olp_\ab^{-1})_{kl} \braket{\ab_l}{\aao_i}\\
                &= \sum_k \ket{\ab_k} (olp_\ab^{-1} olp_{\ab,\aao})_{ki}\\
            \proj_\aao \ket{\ab_i}
                &= \sum_{kl} \ket{\aao_k} (olp_aao^{-1})_{kl} \braket{\aao_l}{\ab_i}\\
                &= \sum_k \ket{aao_k} (olp_\aao^{-1} olp_{\aao,\ab})_{ki}\\
            \proj_\dmo \ket{\ab_i}
                &= \sum_{k} \ket{\dmo_k} \braket{\dmo_k}{\ab_i}\\
                &= \sum_{k} \proj_{\ab} \proj_{aao} \ket{\mo_k}\\
                   \bra{\mo_k} \proj_\aao \proj_\ab \ket{\ab_i}\\
                &= \sum_{k} \ket{\ab_k}
                   (olp_{\ab,\aao} olp_\aao^{-1} olp_{\aao,\mo} olp_{\mo,\aao} olp_\aao^{-1}
                    olp_{\aao,\ab} olp_\ab^{-1} olp_{\ab,\ab})_{ki}\\
            \proj_\mo \ket{\ab_i}
                &= \sum_k \ket{\mo_k} \braket{\mo_k}{\ab_i}\\
            \proj_\dmo \ket{\aao_i}
                &= \sum_{k} \ket{\dmo_k} \braket{\dmo_k}{\aao_i}\\
                &= \sum_{k} \proj_{\ab} \proj_{aao} \ket{\mo_k}
                   \bra{\mo_k} \proj_\aao \proj_\ab \ket{\aao_i}\\
                &= \sum_{k} \ket{\ab_k} (olp_{\ab,\aao} olp_\aao^{-1} olp_{\aao,\mo} olp_{\mo,\aao}
                   olp_\aao^{-1} olp_{\aao,\ab} olp_\ab^{-1} olp_{\ab,\aao})_{ki}\\
            \proj_\mo \ket{\aao_i}
                &= \sum_k \ket{\mo_k} \braket{\mo_k}{\aao_i}\\
                &= \sum_{l} \ket{\ab_l} (C olp_{\mo,\aao})_{li}

        iao type 1:
        p4835
        .. math::

            \ket{iao_i}
                &= (\proj_\omo \proj_\domo + (1-\proj_\omo) (1-\proj_domo)) \proj_\ab \ket{\aao_i}
                &= (1 - \proj_\omo - \proj_\domo + 2 \proj_\omo \proj_\domo) \proj_\ab \ket{\aao_i}

        iao type 2:
        p4841 (in appendix)
        .. math::

            \ket{\iao_i} = ( 1+\proj_\omo-\proj_\domo ) \ket{\aao_i}

        Note that the :math:`\ket{\iao_i}` depends on both :math:`\ab` and :math:`\aao`


        References
        ----------
        .. [1] Knizia, G. Intrinsic Atomic Orbitals: An Unbiased Bridge between
           Quantum Theory and Chemical Concepts. J. Chem. Theory Comput. 2013,
           9, 4834-4843.

        """
        # Get the depolarized occupied mo's
        #  projection of omo onto aao, expressed wrt aao
        coeff_aao_projomo = project(self.olp_aao_aao, self.olp_aao_omo)
        #  projection of aao onto ab, expressed wrt ab
        coeff_ab_projaao = project(self.olp_ab_ab, self.olp_ab_aao)
        #  projecting the omo onto aao, then onto ab
        coeff_ab_domo = coeff_ab_projaao.dot(coeff_aao_projomo)
        #  normalize
        olp_domo_domo = coeff_ab_domo.T.dot(self.olp_ab_ab).dot(coeff_ab_domo)
        coeff_ab_domo *= np.diag(olp_domo_domo) ** (-0.5)
        #  overlap
        olp_domo_domo = coeff_ab_domo.T.dot(self.olp_ab_ab).dot(coeff_ab_domo)
        olp_domo_ab = coeff_ab_domo.T.dot(self.olp_ab_ab)
        if iaotype == 1:
            """ IAO_i = ( O \tilde{O} + (1-O)(1-\tilde{O}) ) proj_ab AAO_i
                      = ( 1 - O - \tilde{O} + 2 O \tilde{O}) proj_ab AAO_i
            where O is the projection onto occupied MO expressed wrt AB
                  \tilde{O} is the projection onto depolarized occupied MO expressed
                               wrt AB
                  proj_ab is the projection onto the atomic basis functions
            """
            # Get projection operators
            def project_onto_omo(coeff_ab_basis):
                """ Projects some basis set onto the occupied molecular orbital space

                Parameter
                ---------
                coeff_ab_basis : np.ndarray(K,M)
                    Transformation from the atomic basis functions to the basis

                Returns
                -------
                coeff_ab_proj : np.ndarray(K,M)
                    Transformation matrix from atomic basis functions to the projection
                    of the basis onto the omo's
                """
                olp_omo_basis = self.olp_omo_ab.dot(coeff_ab_basis)
                return self.coeff_ab_omo.dot(project(self.olp_omo_omo, olp_omo_basis))

            def project_onto_domo(coeff_ab_basis):
                """ Projects some basis set onto the depolarized occupied molecular
                orbital space

                Parameter
                ---------
                coeff_ab_basis : np.ndarray(K,M)
                    Transformation from the atomic basis functions to the basis

                Returns
                -------
                coeff_ab_proj : np.ndarray(K,M)
                    Transformation matrix from atomic basis functions to the projection
                    of the basis onto the omo's
                """
                olp_domo_basis = olp_domo_ab.dot(coeff_ab_basis)
                return coeff_ab_domo.dot(project(olp_domo_domo, olp_domo_basis))

            coeff_basis_iao = np.copy(coeff_ab_projaao)
            coeff_basis_iao -= project_onto_omo(coeff_ab_projaao)
            coeff_basis_iao -= project_onto_domo(coeff_ab_projaao)
            coeff_basis_iao += 2 * project_onto_omo(project_onto_domo(coeff_ab_projaao))
            olp_basis_basis = self.olp_ab_ab
        elif iaotype == 2:
            """ IAO_i = (1 + O - \tilde{O}) AAO_i
                      = AAO_i + (O - \tilde{O}) AAO_i
            where O is the projection onto occupied MO expressed wrt AB
                  \tilde{O} is the projection onto depolarized occupied MO expressed
                               wrt AB

            Note that IAO now depends on both the AAO's and the AB's
            """
            # we can't use the projection operators above, because they operate
            # wrt to the atomic basis functions

            # Project AAO onto OMO's (expressed wrt AB)
            coeff_omo_projaao = project(self.olp_omo_omo, self.olp_omo_aao)
            #  express wrt ab
            coeff_ab_omo_projaao = self.coeff_ab_omo.dot(coeff_omo_projaao)

            # Project AAO onto DOMO's (expressed wrt AB)
            olp_domo_aao = coeff_ab_domo.T.dot(self.olp_ab_aao)
            coeff_domo_projaao = project(olp_domo_domo, olp_domo_aao)
            coeff_ab_domo_projaao = coeff_ab_domo.dot(coeff_domo_projaao)

            # Express IAO wrt NEW basis
            print(
                "WARNING: This transformation matrix is expressed with respect"
                " to BOTH AB and AAO's."
            )
            coeff_basis_iao = np.vstack(
                (coeff_ab_omo_projaao - coeff_ab_domo_projaao, np.identity(self.num_aao))
            )
            olp_basis_basis = np.zeros([self.num_ab + self.num_aao] * 2)
            olp_basis_basis[: self.num_ab] = np.hstack((self.olp_ab_ab, self.olp_ab_aao))
            olp_basis_basis[self.num_ab :] = np.hstack((self.olp_aao_ab, self.olp_aao_aao))
        else:
            raise QuasiTransformationError("Unsupported iaotype, {0}".format(iaotype))
        # orthogonalize
        if is_orth:
            olp_iao_iao = coeff_basis_iao.T.dot(olp_basis_basis).dot(coeff_basis_iao)
            coeff_basis_iao = coeff_basis_iao.dot(orth.power_symmetric(olp_iao_iao, -0.5))
            if not np.allclose(
                np.diag(coeff_basis_iao.T.dot(olp_basis_basis).dot(coeff_basis_iao)), 1, atol=1e-2
            ):
                raise QuasiTransformationError("Cannot orthonormalize IAO's.")
        # normalize
        olp_iao_iao = coeff_basis_iao.T.dot(olp_basis_basis).dot(coeff_basis_iao)
        coeff_basis_iao *= np.diag(olp_iao_iao) ** (-0.5)
        return coeff_basis_iao

    def simple(self, is_orth=False):
        """Construct a basis by projecting the AAO's onto the occupied MO's.

        Returns
        -------
        coeff_ab_simple : np.ndarray
            Transformation from atomic basis to the "simple" basis functions

        Raises
        ------
        QuasiTransformationError
            If the simple projected basis set cannot be orthonormalized.
        """
        coeff_omo_proj = project(self.olp_omo_omo, self.olp_omo_aao)
        # orthogonalize
        if is_orth:
            olp_proj_proj = coeff_omo_proj.T.dot(self.olp_omo_omo).dot(coeff_omo_proj)
            coeff_omo_proj = coeff_omo_proj.dot(orth.power_symmetric(olp_proj_proj, -0.5))
            if not np.allclose(
                np.diag(coeff_omo_proj.T.dot(self.olp_omo_omo).dot(coeff_omo_proj)), 1, atol=1e-2
            ):
                raise QuasiTransformationError("Projected basis set cannot be orthonormalized.")
        # normalize
        olp_proj_proj = coeff_omo_proj.T.dot(self.olp_omo_omo).dot(coeff_omo_proj)
        coeff_omo_proj *= np.diag(olp_proj_proj) ** (-0.5)
        return self.coeff_ab_omo.dot(coeff_omo_proj)


def project(olp_A_A, olp_A_B):
    """Project one basis set (B) onto another basis set (A).

    .. math::

    proj_A \ket{b_i} &=  \sum_{kl} \ket{a_k} (S_{A}^{-1})_{kl} \braket{a_l}{b_i}
    &= \sum_{kl} \ket{a_k} (S_{A}^{-1})_{kl} (S_{A,B})_{li}
    &= \sum_k \ket{a_k} C_{ki}

    where :math:`C_{ki} = \sum_l (S_{A}^{-1})_{kl} (S_{A,B})_{li}`

    Parameters
    ----------
    olp_A_A : np.ndarray(N,N)
        Overlap of the basis set A with basis set A
    olp_A_B : np.ndarray(N,M)
        Overlap of the basis set A with basis set B

    Returns
    -------
    coeff : np.ndarray(N,M)
        Transformation matrix from basis set A to the projection of B onto A
    """
    olp_A_A_inv = orth.power_symmetric(olp_A_A, -1)
    coeff_A_proj = olp_A_A_inv.dot(olp_A_B)
    # Remove zero columns
    coeff_A_proj = coeff_A_proj[:, np.any(coeff_A_proj, axis=0)]
    # Normalize
    olp_proj_proj = coeff_A_proj.T.dot(olp_A_A).dot(coeff_A_proj)
    normalizer = np.diag(olp_proj_proj) ** (-0.5)
    coeff_A_proj *= normalizer
    # Check linear dependence
    rank = np.linalg.matrix_rank(coeff_A_proj)
    if rank < coeff_A_proj.shape[1]:
        print(
            "Warning: There are {0} linearly dependent projections. "
            "The transformation matrix has a shape of {1} and rank of {2}"
            "".format(coeff_A_proj.shape[1] - rank, coeff_A_proj.shape, rank)
        )
    return coeff_A_proj
