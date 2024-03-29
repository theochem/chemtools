"""Tests for orbtools.quasi."""
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path

import numpy as np
from numpy.testing import assert_raises
from chemtools.orbstools.partition import OrbitalPartitionTools
from chemtools.orbstools.quasi import _check_input, make_mmo, project, quambo, quao


def test_project():
    """Test the orbstools.quasi.project."""
    olp_1 = np.identity(10)
    # (trivial) projecting onto same space
    olp_1_2 = np.identity(10)
    assert np.allclose(project(olp_1, olp_1_2), np.identity(10))
    # (trivial) projecting onto subspace
    olp_1_2 = np.eye(10, 6)
    assert np.allclose(project(olp_1, olp_1_2), np.eye(10, 6))
    # (trivial) projecting onto same space + orthogonal complement
    olp_1_2 = np.eye(10, 20)
    assert np.allclose(project(olp_1, olp_1_2), np.eye(10, 10))
    # projecting onto non normalized functions
    olp_1_2 = np.identity(10) * 2
    assert np.allclose(project(olp_1, olp_1_2), np.identity(10))
    # projecting onto linearly dependent basis functions
    olp_1 = np.identity(20)
    olp_1[:10, 10:] = np.identity(10)
    olp_1[10:, :10] = np.identity(10)
    olp_1_2 = np.vstack([np.identity(10)] * 2)
    coeff_1_proj = project(olp_1, olp_1_2)
    assert np.allclose(olp_1.dot(coeff_1_proj), olp_1_2)
    # projecting linearly dependent basis functions
    olp_1 = np.identity(10)
    olp_1_2 = np.hstack([np.identity(10)] * 2)
    assert np.allclose(project(olp_1, olp_1_2), np.hstack([np.identity(10)] * 2))
    # errors
    assert_raises(TypeError, project, olp_1.tolist(), olp_1_2)
    assert_raises(TypeError, project, olp_1.reshape(10, 10, 1), olp_1_2)
    assert_raises(TypeError, project, olp_1.reshape(4, 25), olp_1_2)
    assert_raises(TypeError, project, olp_1, olp_1_2.tolist())
    assert_raises(TypeError, project, olp_1, olp_1_2.reshape(10, 20, 1))
    assert_raises(ValueError, project, olp_1, olp_1_2.T)


def normalize(olp, coeff):
    """Helper function for normalizing orbital coefficients."""
    norm = np.diag(coeff.T.dot(olp).dot(coeff))
    return coeff * norm ** (-0.5)


def test_check_input():
    """Test the orbstools.quasi._check_input."""
    # pylint: disable=R0915
    # olp_ab_ab
    coeff_prim_ab = np.random.rand(10, 10)
    olp_ab_ab = coeff_prim_ab.T.dot(np.identity(10)).dot(coeff_prim_ab)
    norm_ab = np.diag(olp_ab_ab) ** (-0.5)
    olp_ab_ab *= norm_ab[:, None] * norm_ab[None, :]
    # olp_aao_aao
    coeff_prim_ab = np.random.rand(5, 5)
    olp_aao_aao = coeff_prim_ab.T.dot(np.identity(5)).dot(coeff_prim_ab)
    norm_aao = np.diag(olp_aao_aao) ** (-0.5)
    olp_aao_aao *= norm_aao[:, None] * norm_aao[None, :]
    # olp_aao_ab
    olp_aao_ab = np.random.rand(5, 10)
    # coeff_ab_mo
    coeff_ab_mo = olp_ab_ab.dot(np.random.rand(10, 10))
    coeff_ab_mo *= np.diag(coeff_ab_mo.T.dot(olp_ab_ab).dot(coeff_ab_mo)) ** (-0.5)
    # indices_span
    indices_span = np.array([True] * 2 + [False] * 8)

    assert_raises(TypeError, _check_input, coeff_ab_mo=coeff_ab_mo.tolist())
    assert_raises(TypeError, _check_input, coeff_ab_mo=coeff_ab_mo.reshape(10, 10, 1))
    _check_input(coeff_ab_mo=coeff_ab_mo)

    assert_raises(TypeError, _check_input, olp_ab_ab=olp_ab_ab.tolist())
    assert_raises(TypeError, _check_input, olp_ab_ab=olp_ab_ab.reshape(10, 10, 1))
    assert_raises(TypeError, _check_input, olp_ab_ab=olp_ab_ab.reshape(4, 25))
    bad_olp_ab_ab = olp_ab_ab.copy()
    bad_olp_ab_ab *= np.arange(10)[:, None]
    bad_olp_ab_ab *= np.arange(10)[None, :]
    assert_raises(ValueError, _check_input, olp_ab_ab=bad_olp_ab_ab)
    bad_olp_ab_ab = np.random.rand(10, 10)
    bad_olp_ab_ab[np.arange(10), np.arange(10)] = 1
    assert_raises(ValueError, _check_input, olp_ab_ab=bad_olp_ab_ab)
    bad_olp_ab_ab = np.random.rand(10, 10)
    bad_olp_ab_ab += bad_olp_ab_ab.T
    bad_olp_ab_ab[np.arange(10), np.arange(10)] = 1
    assert_raises(ValueError, _check_input, olp_ab_ab=bad_olp_ab_ab)
    _check_input(olp_ab_ab=olp_ab_ab)

    assert_raises(TypeError, _check_input, olp_aao_ab=olp_aao_ab.tolist())
    assert_raises(TypeError, _check_input, olp_aao_ab=olp_aao_ab.reshape(5, 10, 1))
    _check_input(olp_aao_ab=olp_aao_ab)

    assert_raises(TypeError, _check_input, olp_aao_aao=olp_aao_aao.tolist())
    assert_raises(TypeError, _check_input, olp_aao_aao=olp_aao_aao.reshape(5, 5, 1))
    assert_raises(TypeError, _check_input, olp_aao_aao=np.random.rand(4, 5))
    assert_raises(ValueError, _check_input, olp_aao_aao=olp_aao_aao * np.random.rand(5))
    bad_olp_aao_aao = np.random.rand(10, 10)
    bad_olp_aao_aao[np.arange(10), np.arange(10)] = 1
    assert_raises(ValueError, _check_input, olp_aao_aao=bad_olp_aao_aao)
    bad_olp_aao_aao = np.random.rand(10, 10)
    bad_olp_aao_aao += bad_olp_aao_aao.T
    bad_olp_aao_aao[np.arange(10), np.arange(10)] = 1
    assert_raises(ValueError, _check_input, olp_aao_aao=bad_olp_aao_aao)
    _check_input(olp_aao_aao=olp_aao_aao)

    assert_raises(
        ValueError, _check_input, coeff_ab_mo=np.random.rand(9, 10), olp_aao_ab=olp_aao_ab
    )
    _check_input(coeff_ab_mo=coeff_ab_mo, olp_aao_ab=olp_aao_ab)

    assert_raises(ValueError, _check_input, olp_ab_ab=olp_ab_ab, olp_aao_ab=np.random.rand(5, 9))
    _check_input(olp_ab_ab=olp_ab_ab, olp_aao_ab=olp_aao_ab)

    assert_raises(ValueError, _check_input, coeff_ab_mo=np.random.rand(9, 10), olp_ab_ab=olp_ab_ab)
    _check_input(coeff_ab_mo=coeff_ab_mo, olp_ab_ab=olp_ab_ab)

    assert_raises(
        ValueError, _check_input, olp_aao_ab=np.random.rand(4, 10), olp_aao_aao=olp_aao_aao
    )
    _check_input(olp_aao_ab=olp_aao_ab, olp_aao_aao=olp_aao_aao)

    assert_raises(ValueError, _check_input, coeff_ab_mo=np.random.rand(10, 10), olp_ab_ab=olp_ab_ab)
    _check_input(coeff_ab_mo=coeff_ab_mo, olp_ab_ab=olp_ab_ab)

    assert_raises(TypeError, _check_input, indices_span=indices_span.tolist())
    assert_raises(TypeError, _check_input, indices_span=indices_span.reshape(10, 1))
    assert_raises(TypeError, _check_input, indices_span=indices_span.astype(int))
    _check_input(indices_span=indices_span)

    bad_indices_span = np.random.rand(9) < 0.5
    assert_raises(ValueError, _check_input, indices_span=bad_indices_span, coeff_ab_mo=coeff_ab_mo)
    _check_input(indices_span=indices_span, coeff_ab_mo=coeff_ab_mo)


def test_make_mmo():
    """Test orbstools.quasi.make_mmo."""
    # olp_aao_ab
    olp_aao_ab = np.random.rand(5, 10)
    # coeff_ab_mo
    olp_ab_ab, _, _ = np.linalg.svd(np.random.rand(10, 10))
    olp_ab_ab = (olp_ab_ab * np.random.rand(10)).dot(olp_ab_ab.T)
    norm_ab = np.diag(olp_ab_ab) ** (-0.5)
    olp_ab_ab *= norm_ab[:, None] * norm_ab[None, :]
    coeff_ab_mo = olp_ab_ab.dot(np.random.rand(10, 10))
    coeff_ab_mo *= np.diag(coeff_ab_mo.T.dot(olp_ab_ab).dot(coeff_ab_mo)) ** (-0.5)
    # indices_span
    indices_span = np.array([True] * 5 + [False] * 5)

    assert_raises(TypeError, make_mmo, olp_aao_ab, coeff_ab_mo, indices_span, dim_mmo=8.0)
    assert_raises(ValueError, make_mmo, olp_aao_ab, coeff_ab_mo, indices_span, dim_mmo=11)
    assert_raises(ValueError, make_mmo, olp_aao_ab, coeff_ab_mo, indices_span, dim_mmo=4)

    coeff_ab_mmo = make_mmo(olp_aao_ab, coeff_ab_mo, indices_span, dim_mmo=6)
    # check that occupied mo's are spanned exactly
    coeff_mmo_newmo = project(
        coeff_ab_mmo.T.dot(olp_ab_ab).dot(coeff_ab_mmo),
        coeff_ab_mmo.T.dot(olp_ab_ab).dot(coeff_ab_mo),
    )
    olp_mo_newmo = coeff_ab_mo.T.dot(olp_ab_ab).dot(coeff_ab_mmo).dot(coeff_mmo_newmo)
    olp_mo_mo = coeff_ab_mo.T.dot(olp_ab_ab).dot(coeff_ab_mo)
    assert np.allclose(np.diag(olp_mo_newmo)[:5], 1)
    assert np.allclose(olp_mo_newmo[:5, :5], olp_mo_mo[:5, :5])
    # check that the virtual mo's are not spanned exactly
    assert not np.allclose(np.diag(olp_mo_newmo)[5:], 1)
    assert not np.allclose(olp_mo_newmo, olp_mo_mo)


def test_make_mmo_old_code():
    """Test orbstools.quasi.make_mmo against output from old code.

    Old code can be found in https://github.com/QuantumElephant/dumbo/tree/master/quasibasis.

    """
    # compare against reference generated using old code
    with path("chemtools.data", "naclo4_olp_aao_ab.npy") as fname:
        olp_aao_ab = np.load(str(fname))
    with path("chemtools.data", "naclo4_coeff_ab_mo.npy") as fname:
        coeff_ab_mo = np.load(str(fname))
    with path("chemtools.data", "naclo4_occupations.npy") as fname:
        occupations = np.load(str(fname))
    indices_span = occupations > 0

    with path("chemtools.data", "naclo4_coeff_ab_mmo.npy") as fname:
        coeff_ab_mmo = np.load(str(fname))
    assert np.allclose(coeff_ab_mmo, make_mmo(olp_aao_ab, coeff_ab_mo, indices_span))


def test_quambo_old_code():
    """Test orbstools.quasi.quambo against output from old code.

    Old code can be found in https://github.com/QuantumElephant/dumbo/tree/master/quasibasis.

    """
    with path("chemtools.data", "naclo4_coeff_ab_mo.npy") as fname:
        coeff_ab_mo = np.load(str(fname))
    with path("chemtools.data", "naclo4_olp_ab_ab.npy") as fname:
        olp_ab_ab = np.load(str(fname))
    with path("chemtools.data", "naclo4_olp_aao_ab.npy") as fname:
        olp_aao_ab = np.load(str(fname))
    with path("chemtools.data", "naclo4_occupations.npy") as fname:
        occupations = np.load(str(fname))
    indices_span = occupations > 0

    with path("chemtools.data", "naclo4_coeff_ab_quambo.npy") as fname:
        coeff_ab_quambo = np.load(str(fname))
    with path("chemtools.data", "naclo4_olp_aao_ab.npy") as fname:
        olp_aao_ab = np.load(str(fname))
    assert np.allclose(coeff_ab_quambo, quambo(olp_ab_ab, olp_aao_ab, coeff_ab_mo, indices_span))


def test_quao_old_code():
    """Test orbstools.quasi.quao against output from old code.

    Old code can be found in https://github.com/QuantumElephant/dumbo/tree/master/quasibasis.

    """
    with path("chemtools.data", "naclo4_coeff_ab_mo.npy") as fname:
        coeff_ab_mo = np.load(str(fname))
    with path("chemtools.data", "naclo4_olp_ab_ab.npy") as fname:
        olp_ab_ab = np.load(str(fname))
    with path("chemtools.data", "naclo4_olp_aao_ab.npy") as fname:
        olp_aao_ab = np.load(str(fname))
    with path("chemtools.data", "naclo4_olp_aao_aao.npy") as fname:
        olp_aao_aao = np.load(str(fname))
    with path("chemtools.data", "naclo4_occupations.npy") as fname:
        occupations = np.load(str(fname))
    indices_span = occupations > 0

    with path("chemtools.data", "naclo4_coeff_ab_quao.npy") as fname:
        coeff_ab_quao = np.load(str(fname))
    with path("chemtools.data", "naclo4_olp_aao_ab.npy") as fname:
        olp_aao_ab = np.load(str(fname))
    assert np.allclose(
        coeff_ab_quao, quao(olp_ab_ab, olp_aao_ab, olp_aao_aao, coeff_ab_mo, indices_span)
    )


def test_quambo():
    """Test orbstools.quasi.quambo against literature value for Mulliken populations.

    References
    ----------
    .. [1] Janowski, T. Near equivalence of intrinsic atomic orbitals and quasiatomic orbitals.
        JCTC, 2014, 10, 3085-3091.

    """
    with path("chemtools.data", "naclo4_coeff_ab_mo.npy") as fname:
        coeff_ab_mo = np.load(str(fname))
    with path("chemtools.data", "naclo4_olp_ab_ab.npy") as fname:
        olp_ab_ab = np.load(str(fname))
    with path("chemtools.data", "naclo4_olp_aao_ab.npy") as fname:
        olp_aao_ab = np.load(str(fname))
    with path("chemtools.data", "naclo4_occupations.npy") as fname:
        occupations = np.load(str(fname))
    indices_span = occupations > 0
    with path("chemtools.data", "naclo4_qab_atom_indices.npy") as fname:
        ab_atom_indices = np.load(str(fname))

    olp_ab_omo = olp_ab_ab.dot(coeff_ab_mo[:, indices_span])

    coeff_ab_quambo = quambo(olp_ab_ab, olp_aao_ab, coeff_ab_mo, indices_span)
    olp_quambo_quambo = coeff_ab_quambo.T.dot(olp_ab_ab).dot(coeff_ab_quambo)
    olp_quambo_omo = coeff_ab_quambo.T.dot(olp_ab_omo)
    coeff_quambo_omo = project(olp_quambo_quambo, olp_quambo_omo)
    pop = OrbitalPartitionTools(
        coeff_quambo_omo, occupations[indices_span], olp_quambo_quambo, 6, ab_atom_indices
    ).mulliken_populations()
    partial_pop = np.array([11, 17, 8, 8, 8, 8]) - pop
    assert np.allclose(
        partial_pop, np.array([0.968, 2.454, -0.807, -0.904, -0.904, -0.807]), atol=1e-3
    )


def test_quao():
    """Test orbstools.quasi.quambo against literature value for Mulliken populations.

    References
    ----------
    .. [1] Janowski, T. Near equivalence of intrinsic atomic orbitals and quasiatomic orbitals.
        JCTC, 2014, 10, 3085-3091.

    """
    with path("chemtools.data", "naclo4_coeff_ab_mo.npy") as fname:
        coeff_ab_mo = np.load(str(fname))
    with path("chemtools.data", "naclo4_olp_ab_ab.npy") as fname:
        olp_ab_ab = np.load(str(fname))
    with path("chemtools.data", "naclo4_olp_aao_ab.npy") as fname:
        olp_aao_ab = np.load(str(fname))
    with path("chemtools.data", "naclo4_olp_aao_aao.npy") as fname:
        olp_aao_aao = np.load(str(fname))
    with path("chemtools.data", "naclo4_occupations.npy") as fname:
        occupations = np.load(str(fname))
    indices_span = occupations > 0
    with path("chemtools.data", "naclo4_qab_atom_indices.npy") as fname:
        ab_atom_indices = np.load(str(fname))

    olp_ab_omo = olp_ab_ab.dot(coeff_ab_mo[:, indices_span])

    coeff_ab_quao = quao(olp_ab_ab, olp_aao_ab, olp_aao_aao, coeff_ab_mo, indices_span)
    olp_quao_quao = coeff_ab_quao.T.dot(olp_ab_ab).dot(coeff_ab_quao)
    olp_quao_omo = coeff_ab_quao.T.dot(olp_ab_omo)
    coeff_quao_omo = project(olp_quao_quao, olp_quao_omo)
    pop = OrbitalPartitionTools(
        coeff_quao_omo, occupations[indices_span], olp_quao_quao, 6, ab_atom_indices
    ).mulliken_populations()
    partial_pop = np.array([11, 17, 8, 8, 8, 8]) - pop
    assert np.allclose(
        partial_pop, np.array([0.967, 2.498, -0.819, -0.914, -0.914, -0.819]), atol=1e-3
    )
