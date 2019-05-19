"""Module for making Quasiatomic orbitals."""
import numpy as np
from chemtools.orbstools import orthogonalization as orth


def _check_input(
    olp_ab_ab=None, olp_aao_ab=None, olp_aao_aao=None, coeff_ab_mo=None, indices_span=None
):
    """Check the inputs.

    Parameters
    ----------
    olp_ab_ab : np.ndarray(K, K)
        Overlaps of the atomic basis functions.
        :math:`K` is the number of atomic basis functions.
    olp_aao_ab : np.ndarray(L, K)
        Overlaps of the atomic basis functions with the reference basis functions (aao).
        Rows correspond to the reference basis functions. :math:`L` is the number of reference
        basis functions.
        Column correspond to the atomic basis functions. :math:`K` is the number of atomic basis
        functions.
    olp_aao_aao : np.ndarray(L, L)
        Overlaps of the atomic basis functions with the reference basis functions (aao).
        :math:`L` is the number of reference basis functions.
    coeff_ab_mo : np.ndarray(K, M)
        Transformation matrix from the atomic basis functions to molecular orbitals. The matrix is
        applied onto the right side.
        Rows correspond to the atomic basis functions. :math:`K` is the number of atomic basis
        functions.
        Columns correspond to the molecular basis functions. :math:`L` is the number of molecular
        basis functions.
    indices_span : np.ndarray(M)
        Molecular orbitals that will be spanned exactly by the quasi basis functions.
        Each entry is a boolean, where molecular orbitals that are exactly described have value
        `True`.

    Raises
    ------
    TypeError
        If `coeff_ab_mo` is two-dimensional numpy array.
        If `olp_ab_ab` is not a two-dimensional square numpy array.
        If `olp_aao_ab` is not a two-dimensional numpy array.
        If `olp_aao_aao` is not a two-dimensional square numpy array.
        If `indices_span` is not a one-dimensional numpy array of dtype bool.
    ValueError
        If `olp_ab_ab` is not normalized, i.e. diagonal that is not equal to 1.
        If `olp_ab_ab` is not symmetric.
        If `olp_ab_ab` is not positive-semidefinite.
        If `olp_aao_aao` is not normalized, i.e. diagonal that is not equal to 1.
        If `olp_aao_aao` is not symmetric.
        If `olp_aao_aao` is not positive-semidefinite.
        If number of atomic basis functions is not consistent between `coeff_ab_mo`, `olp_ab_ab`,
        and `olp_aao_ab`.
        If number of reference basis functions (aao) is not consistent for `olp_aao_ab` and
        `olp_aao_aao`.
        If molecular orbitals are not normalized.

    """
    # pylint: disable=R0912
    if coeff_ab_mo is not None:
        if not (isinstance(coeff_ab_mo, np.ndarray) and coeff_ab_mo.ndim == 2):
            raise TypeError("Given coefficient matrix is not a two-dimensional numpy array.")

    if olp_ab_ab is not None:
        if not (
            isinstance(olp_ab_ab, np.ndarray)
            and olp_ab_ab.ndim == 2
            and olp_ab_ab.shape[0] == olp_ab_ab.shape[1]
        ):
            raise TypeError(
                "Given overlap matrix for atomic basis is not a two-dimensional square numpy array."
            )
        if not np.allclose(np.diag(olp_ab_ab), np.ones(olp_ab_ab.shape[0])):
            raise ValueError("Given overlap matrix for atomic basis is not normalized.")
        if not np.allclose(olp_ab_ab, olp_ab_ab.T):
            raise ValueError("Given overlap matrix for atomic basis is not symmetric.")
        if not np.all(orth.eigh(olp_ab_ab)[0] >= 0):
            raise ValueError("Given overlap matrix for atomic basis is not positive semidefinite.")

    if olp_aao_ab is not None:
        if not (isinstance(olp_aao_ab, np.ndarray) and olp_aao_ab.ndim == 2):
            raise TypeError(
                "Given overlap matrix for atomic basis and AAO is not a two-dimensional numpy "
                "array."
            )

    if olp_aao_aao is not None:
        if not (
            isinstance(olp_aao_aao, np.ndarray)
            and olp_aao_aao.ndim == 2
            and olp_aao_aao.shape[0] == olp_aao_aao.shape[1]
        ):
            raise TypeError(
                "Given overlap matrix for AAO is not a two dimensional square numpy array."
            )
        if not np.allclose(np.diag(olp_aao_aao), np.ones(olp_aao_aao.shape[0])):
            raise ValueError("Given overlap matrix for AAO is not normalized.")
        if not np.allclose(olp_aao_aao, olp_aao_aao.T):
            raise ValueError("Given overlap matrix for AAO is not symmetric.")
        if not np.all(orth.eigh(olp_aao_aao)[0] >= 0):
            raise ValueError("Given overlap matrix for AAO is not positive semidefinite.")

    if (
        coeff_ab_mo is not None
        and olp_aao_ab is not None
        and coeff_ab_mo.shape[0] != olp_aao_ab.shape[1]
    ):
        raise ValueError(
            "Number of rows of the `coeff_ab_mo` is not equal to the number of columns of "
            "`olp_aao_ab`."
        )
    if (
        olp_ab_ab is not None
        and olp_aao_ab is not None
        and olp_ab_ab.shape[0] != olp_aao_ab.shape[1]
    ):
        raise ValueError(
            "Number of rows of the `olp_ab_ab` is not equal to the number of columns of "
            "`olp_aao_ab`."
        )
    if (
        coeff_ab_mo is not None
        and olp_ab_ab is not None
        and coeff_ab_mo.shape[0] != olp_ab_ab.shape[0]
    ):
        raise ValueError(
            "Number of rows of the `coeff_ab_mo` is not equal to the number of rows/columns of "
            "`olp_ab_ab`."
        )

    if (
        olp_aao_ab is not None
        and olp_aao_aao is not None
        and olp_aao_ab.shape[0] != olp_aao_aao.shape[0]
    ):
        raise ValueError(
            "Number of rows of the `olp_aao_ab` is not equal to the number of rows/columns of "
            "`olp_ab_ab`."
        )

    if coeff_ab_mo is not None and olp_ab_ab is not None:
        olp_mo_mo = coeff_ab_mo.T.dot(olp_ab_ab).dot(coeff_ab_mo)
        if not np.allclose(np.diag(olp_mo_mo), np.ones(olp_mo_mo.shape[0])):
            raise ValueError(
                "The overlap of the molecular orbitals, calculated from `coeff_ab_mo` and "
                "`olp_ab_ab` is not normalized."
            )

    if indices_span is not None:
        if not (
            isinstance(indices_span, np.ndarray)
            and indices_span.ndim == 1
            and indices_span.dtype == bool
        ):
            raise TypeError(
                "`indices_span` must be given as a one-dimensional numpy array of dtype bool."
            )
        if coeff_ab_mo is not None and indices_span.size != coeff_ab_mo.shape[1]:
            raise ValueError(
                "`indices_span` must have as many entries as there are columns in the "
                "`coeff_ab_mo` (i.e. number of molecular orbitals)."
            )


def project(olp_one_one, olp_one_two):
    r"""Project one basis set onto another basis set.

    .. math::

        \mathrm{proj}_A \ket{b_i} &=  \sum_{kl} \ket{a_k} (S_{A}^{-1})_{kl} \braket{a_l | b_i}\\
        &= \sum_{kl} \ket{a_k} (S_{A}^{-1})_{kl} (S_{A,B})_{li}\\
        &= \sum_k \ket{a_k} C_{ki}

    where :math:`(S_A)_{kl} = \braket{a_k | a_l}`, :math:`(S_{A,B})_{li} = \braket{a_l | b_i}`, and
    :math:`C_{ki} = \sum_l (S_{A}^{-1})_{kl} (S_{A,B})_{li}`.

    Parameters
    ----------
    olp_one_one : np.ndarray(N, N)
        Overlap of the basis functions in set 1 with basis functions from set 1.
    olp_one_two : np.ndarray(N, M)
        Overlap of the basis functions in set 1 with basis functions from set 2.

    Returns
    -------
    coeff : np.ndarray(N, M)
        Transformation matrix from basis functions in set 1 to the projection of baiss set 2 onto
        basis set 1.

    """
    if not (
        isinstance(olp_one_one, np.ndarray)
        and olp_one_one.ndim == 2
        and olp_one_one.shape[0] == olp_one_one.shape[1]
    ):
        raise TypeError("`olp_one_one` must be a two-dimensional square numpy array.")
    if not (isinstance(olp_one_two, np.ndarray) and olp_one_two.ndim == 2):
        raise TypeError("`olp_one_two` must be a two-dimensional numpy array.")
    if olp_one_one.shape[0] != olp_one_two.shape[0]:
        raise ValueError(
            "Number of rows/columns of `olp_one_one` must be equal to the number of rows in "
            "`olp_one_two`."
        )
    olp_one_one_inv = orth.power_symmetric(olp_one_one, -1)
    coeff_one_proj = olp_one_one_inv.dot(olp_one_two)
    # Remove zero columns
    coeff_one_proj = coeff_one_proj[:, np.any(coeff_one_proj, axis=0)]
    # Normalize
    olp_proj_proj = coeff_one_proj.T.dot(olp_one_one).dot(coeff_one_proj)
    normalizer = np.diag(olp_proj_proj) ** (-0.5)
    coeff_one_proj *= normalizer
    # Check linear dependence
    rank = np.linalg.matrix_rank(coeff_one_proj)
    if rank < coeff_one_proj.shape[1]:
        print(
            "Warning: There are {0} linearly dependent projections. The transformation matrix has a"
            " shape of {1} and rank of {2}".format(
                coeff_one_proj.shape[1] - rank, coeff_one_proj.shape, rank
            )
        )
    return coeff_one_proj


def make_mmo(olp_aao_ab, coeff_ab_mo, indices_span, dim_mmo=None):
    r"""Return transformation matrix from atomic basis functions to minimal molecular orbitals.

    Parameters
    ----------
    olp_aao_ab : np.ndarray(M, N)
        Overlap between reference basis functions (rows) and atomic basis functions (columns).
    coeff_ab_mo : np.ndarray(K, N)
        Transformation matrix from atomic basis functions (rows) to molecular orbitals (columns).
        Transformation is applied to the right sid.
    indices_span : np.ndarray(N)
        Boolean indices for the molecular orbitals that will be spanned by the generated MMO's.
    dim_mmo : {int, None}
        Total dimension of the MMO space.
        Default is the dimension of the reference basis function space.

    Returns
    -------
    coeff_ab_mmo : np.ndarray
        Transformation matrix from atomic basis functions to mMO's.

    Raises
    ------
    TypeError
        If `dim_mmo` is not an integer (or None).
    ValueError
        If the dimension of the MMO space is larger than the number of molecular orbitals.
        If the dimension of the MMO space is smaller than the space that needs to be spanned.

    References
    ----------
    .. [1] Lu. W.C.; Wang, C.Z.; Schmidt, W.; Bytautas, L.;Ho K.M.; Ruedenberg, K. Ruedenberg.
        Molecule intrinsic minimal basis sets. I. Exact resolution of ab initio optimized molecular
        orbitals in terms of deformed atomic minimal.
    .. [2] West, A.C.; Schmidt, M.W. A comprehensive analysis of molecule-intrinsic quasiatomic,
        bonding, and correlating orbitals. I. Hartree-Fock wave functions. J. Chem. Phys. 2013, 139,
        234107.

    """
    _check_input(coeff_ab_mo=coeff_ab_mo, olp_aao_ab=olp_aao_ab, indices_span=indices_span)
    olp_aao_mo = olp_aao_ab.dot(coeff_ab_mo)

    num_aao, num_mo = olp_aao_mo.shape
    if dim_mmo is None:
        dim_mmo = num_aao
    if not isinstance(dim_mmo, int):
        raise TypeError("Dimension of MMO space must be an integer (or None).")
    if dim_mmo > num_mo:
        raise ValueError(
            "Dimension of MMO space, {0}, is larger than the number of molecular orbitals, {1}."
            "".format(dim_mmo, num_mo)
        )
    if dim_mmo < np.sum(indices_span):
        raise ValueError(
            "Dimension of MMO space, {0}, is smaller than the space you want to span, {1}."
            "".format(dim_mmo, np.sum(indices_span))
        )

    # Set local variables
    dim_span = np.sum(indices_span)
    num_to_add = dim_mmo - dim_span

    # Create occupied MMO
    coeff_occmo_occmmo = np.identity(dim_span)

    # Create virtual MMO
    #  find overlap between aao and virtuals
    olp_aao_virmo = olp_aao_mo[:, ~indices_span]
    #  from the right singular vector of olp_aao_virmo
    coeff_virmo_virmmo = orth.svd(olp_aao_virmo)[2].T
    #  select vectors with largest (num_to_add) singular values
    coeff_virmo_virmmo = coeff_virmo_virmmo[:, :num_to_add]

    # Combine the occupied and virtual MMO's
    coeff_mo_occmmo = np.zeros((indices_span.size, coeff_occmo_occmmo.shape[1]))
    coeff_mo_occmmo[indices_span] = coeff_occmo_occmmo
    coeff_mo_virmmo = np.zeros((indices_span.size, coeff_virmo_virmmo.shape[1]))
    coeff_mo_virmmo[~indices_span] = coeff_virmo_virmmo
    coeff_mo_mmo = np.hstack((coeff_mo_occmmo, coeff_mo_virmmo))

    # Express MMO wrt atomic basis functions
    return coeff_ab_mo.dot(coeff_mo_mmo)


def quambo(olp_ab_ab, olp_aao_ab, coeff_ab_mo, indices_span, dim=None):
    r"""Return transformation matrix from atomic basis functions to QUAMBO's.

    Parameters
    ----------
    olp_ab_ab : np.ndarray(K, K)
        Overlaps of the atomic basis functions.
        :math:`K` is the number of atomic basis functions.
    olp_aao_ab : np.ndarray(L, K)
        Overlaps of the reference basis functions (aao) with the atomic basis functions.
        Rows correspond to the atomic basis functions. :math:`K` is the number of atomic basis
        functions.
        Columns correspond to the reference basis functions. :math:`L` is the number of reference
        basis functions.
    coeff_ab_mo : np.ndarray(K, M)
        Transformation matrix from the atomic basis functions to molecular orbitals. The matrix is
        applied onto the right side.
        Rows correspond to the atomic basis functions. :math:`K` is the number of atomic basis
        functions.
        Columns correspond to the molecular basis functions. :math:`L` is the number of molecular
        basis functions.
    indices_span : np.ndarray(M)
        Molecular orbitals that will be spanned exactly by the QUAMBO's.
        Each entry is a boolean, where molecular orbitals that are exactly described have value
        `True`.
    dim : {int, None}
        Number of QUAMBO basis functions.
        Default is the number of reference basis functions.

    Returns
    -------
    coeff_ab_quambo : np.ndarray
        Transformation matrix from atomic basis functions to QUAMBO's.

    Notes
    -----
    .. math::

        \ket{\mathrm{QUAMBO}_k}
        &= \mathrm{proj}_{\mathrm{mMO}^\mathrm{QUAMBO}} \ket{\mathrm{AAO}_k}\\
        &= \sum_i \ket{\mathrm{mMO}^\mathrm{QUAMBO}_i}
            (olp_{{\mathrm{mMO}^\mathrm{QUAMBO}}}^{-1})_{ij}
            \braket{ \mathrm{mMO}^\mathrm{QUAMBO}_j | \mathrm{AAO}_k }\\
        &= \sum_j \ket{AB_j} (C^{\mathrm{AB},\mathrm{mMO}^\mathrm{QUAMBO}}
            olp_{{\mathrm{mMO}^\mathrm{QUAMBO}}}^{-1}
            olp_{\mathrm{mMO}^\mathrm{QUAMBO}, \mathrm{AAO}})_{ji}\\

    References
    ----------
    .. [1] Lu. W.C.; Wang, C.Z.; Schmidt, W.; Bytautas, L.;Ho K.M.; Ruedenberg, K. Ruedenberg.
        Molecule intrinsic minimal basis sets. I. Exact resolution of ab initio optimized molecular
        orbitals in terms of deformed atomic minimal.

    """
    _check_input(
        olp_ab_ab=olp_ab_ab,
        olp_aao_ab=olp_aao_ab,
        coeff_ab_mo=coeff_ab_mo,
        indices_span=indices_span,
    )
    # Find MMO for QUAMBOs
    coeff_ab_mmo = make_mmo(olp_aao_ab, coeff_ab_mo, indices_span, dim_mmo=dim)
    # Get transformation
    olp_mmo_mmo = coeff_ab_mmo.T.dot(olp_ab_ab).dot(coeff_ab_mmo)
    olp_mmo_aao = (olp_aao_ab.dot(coeff_ab_mmo)).T
    coeff_mmo_proj = project(olp_mmo_mmo, olp_mmo_aao)
    # Normalize
    olp_proj_proj = coeff_mmo_proj.T.dot(olp_mmo_mmo).dot(coeff_mmo_proj)
    coeff_mmo_proj *= np.diag(olp_proj_proj) ** (-0.5)

    return coeff_ab_mmo.dot(coeff_mmo_proj)


def quao(olp_ab_ab, olp_aao_ab, olp_aao_aao, coeff_ab_mo, indices_span, dim=None):
    r"""Return transformation matrix from atomic basis functions to QUAO's.

    Parameters
    ----------
    olp_ab_ab : np.ndarray(K, K)
        Overlaps of the atomic basis functions.
        :math:`K` is the number of atomic basis functions.
    olp_ab_aao : np.ndarray(K, L)
        Overlaps of the atomic basis functions with the reference basis functions (aao).
        Rows correspond to the atomic basis functions. :math:`K` is the number of atomic basis
        functions.
        Columns correspond to the reference basis functions. :math:`L` is the number of reference
        basis functions.
    olp_aao_aao : np.ndarray(L, L)
        Overlaps of the atomic basis functions with the reference basis functions (aao).
        Rows correspond to the atomic basis functions. :math:`K` is the number of atomic basis
        functions.
        Columns correspond to the reference basis functions. :math:`L` is the number of reference
        basis functions.
    coeff_ab_mo : np.ndarray(K, M)
        Transformation matrix from the atomic basis functions to molecular orbitals. The matrix is
        applied onto the right side.
        Rows correspond to the atomic basis functions. :math:`K` is the number of atomic basis
        functions.
        Columns correspond to the molecular basis functions. :math:`L` is the number of molecular
        basis functions.
    indices_span : np.ndarray(M)
        Molecular orbitals that will be spanned exactly by the QUAMBO's.
        Each entry is a boolean, where molecular orbitals that are exactly described have value
        `True`.
    dim : {int, None}
        Number of QUAMBO basis functions.
        Default is the number of reference basis functions.

    Returns
    -------
    coeff_ab_quao
        Transformation matrix from atomic basis functions to QUAO's

    Notes
    -----
    .. math::

        \ket{\mathrm{QUAO}_k} &= \mathrm{proj}_{\mathrm{mMO}^\mathrm{QUAO}}\\
        &= \sum_i \ket{ \mathrm{mMO}^\mathrm{QUAO}_i }
            (olp_{{\mathrm{mMO}^\mathrm{QUAO}}}^{-1})_{ij}
            \braket{ \mathrm{mMO}^\mathrm{QUAO}_j | \mathrm{AAO}_k }\\
        &= \sum_j \ket{AB_j} (C^{\mathrm{AB},\mathrm{mMO}^\mathrm{QUAO}}
            olp_{{\mathrm{mMO}^\mathrm{QUAO}}}^{-1}
            olp_{\mathrm{mMO}^\mathrm{QUAO}, \mathrm{AAO}})_{ji}\\

    References
    ----------
    .. [1] West, A.C.; Schmidt, M.W. Gordon, M.S; Ruedenberg, K. A comprehensive analysis of
        molecule-intrinsic quasiatomic, bonding, and correlating orbitals. I. Hartree-Fock wave
        functions. J. Chem. Phys. 2013, 139, 234107.

    """
    _check_input(
        olp_ab_ab=olp_ab_ab,
        olp_aao_ab=olp_aao_ab,
        olp_aao_aao=olp_aao_aao,
        coeff_ab_mo=coeff_ab_mo,
        indices_span=indices_span,
    )
    # Orthogonalize AAOs
    olp_oaao_ab = orth.power_symmetric(olp_aao_aao, -0.5).dot(olp_aao_ab)

    # Get MMOs using the orthogonalized AAOs (MMO for QUAOs)
    coeff_ab_mmo = make_mmo(olp_oaao_ab, coeff_ab_mo, indices_span, dim_mmo=dim)

    # Find transformation for QUAOs
    olp_mmo_mmo = coeff_ab_mmo.T.dot(olp_ab_ab).dot(coeff_ab_mmo)
    olp_mmo_aao = (olp_aao_ab.dot(coeff_ab_mmo)).T
    coeff_mmo_proj = project(olp_mmo_mmo, olp_mmo_aao)

    # Normalize
    olp_proj_proj = coeff_mmo_proj.T.dot(olp_mmo_mmo).dot(coeff_mmo_proj)
    coeff_mmo_proj *= np.diag(olp_proj_proj) ** (-0.5)

    return coeff_ab_mmo.dot(coeff_mmo_proj)
