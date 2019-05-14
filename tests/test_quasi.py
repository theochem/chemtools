"""Tests for orbtools.quasi."""
import os

import numpy as np
from orbtools.mulliken import mulliken_populations
from orbtools.quasi import _check_input, make_mmo, project, quambo, quao
import pytest


def test_project():
    """Test the orbtools.quasi.project."""
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
    with pytest.raises(TypeError):
        project(olp_1.tolist(), olp_1_2)
    with pytest.raises(TypeError):
        project(olp_1.reshape(10, 10, 1), olp_1_2)
    with pytest.raises(TypeError):
        project(olp_1.reshape(4, 25), olp_1_2)
    with pytest.raises(TypeError):
        project(olp_1, olp_1_2.tolist())
    with pytest.raises(TypeError):
        project(olp_1, olp_1_2.reshape(10, 20, 1))
    with pytest.raises(ValueError):
        project(olp_1, olp_1_2.T)


def normalize(olp, coeff):
    norm = np.diag(coeff.T.dot(olp).dot(coeff))
    return coeff * norm ** (-0.5)


def test_check_input():
    """Test the orbtools.quasi._check_input."""
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

    with pytest.raises(TypeError):
        _check_input(coeff_ab_mo=coeff_ab_mo.tolist())
    with pytest.raises(TypeError):
        _check_input(coeff_ab_mo=coeff_ab_mo.reshape(10, 10, 1))
    _check_input(coeff_ab_mo=coeff_ab_mo)

    with pytest.raises(TypeError):
        _check_input(olp_ab_ab=olp_ab_ab.tolist())
    with pytest.raises(TypeError):
        _check_input(olp_ab_ab=olp_ab_ab.reshape(10, 10, 1))
    with pytest.raises(TypeError):
        _check_input(olp_ab_ab=olp_ab_ab.reshape(4, 25))
    with pytest.raises(ValueError):
        _check_input(olp_ab_ab=olp_ab_ab * np.random.rand(10))
    with pytest.raises(ValueError):
        _check_input(olp_ab_ab=np.random.rand(10, 10))
    with pytest.raises(ValueError):
        bad_olp_ab_ab = np.random.rand(10, 10)
        _check_input(olp_ab_ab=bad_olp_ab_ab + bad_olp_ab_ab.T)
    _check_input(olp_ab_ab=olp_ab_ab)

    with pytest.raises(TypeError):
        _check_input(olp_aao_ab=olp_aao_ab.tolist())
    with pytest.raises(TypeError):
        _check_input(olp_aao_ab=olp_aao_ab.reshape(5, 10, 1))
    _check_input(olp_aao_ab=olp_aao_ab)

    with pytest.raises(TypeError):
        _check_input(olp_aao_aao=olp_aao_aao.tolist())
    with pytest.raises(TypeError):
        _check_input(olp_aao_aao=olp_aao_aao.reshape(5, 5, 1))
    with pytest.raises(TypeError):
        _check_input(olp_aao_aao=np.random.rand(4, 5))
    with pytest.raises(ValueError):
        _check_input(olp_aao_aao=olp_aao_aao * np.random.rand(5))
    with pytest.raises(ValueError):
        _check_input(olp_aao_aao=np.random.rand(5, 5))
    with pytest.raises(ValueError):
        bad_olp_aao_aao = np.random.rand(5, 5)
        _check_input(olp_aao_aao=bad_olp_aao_aao + bad_olp_aao_aao.T)
    _check_input(olp_aao_aao=olp_aao_aao)

    with pytest.raises(ValueError):
        _check_input(coeff_ab_mo=np.random.rand(9, 10), olp_aao_ab=olp_aao_ab)
    _check_input(coeff_ab_mo=coeff_ab_mo, olp_aao_ab=olp_aao_ab)

    with pytest.raises(ValueError):
        _check_input(olp_ab_ab=olp_ab_ab, olp_aao_ab=np.random.rand(5, 9))
    _check_input(olp_ab_ab=olp_ab_ab, olp_aao_ab=olp_aao_ab)

    with pytest.raises(ValueError):
        _check_input(coeff_ab_mo=np.random.rand(9, 10), olp_ab_ab=olp_ab_ab)
    _check_input(coeff_ab_mo=coeff_ab_mo, olp_ab_ab=olp_ab_ab)

    with pytest.raises(ValueError):
        _check_input(olp_aao_ab=np.random.rand(4, 10), olp_aao_aao=olp_aao_aao)
    _check_input(olp_aao_ab=olp_aao_ab, olp_aao_aao=olp_aao_aao)

    with pytest.raises(ValueError):
        _check_input(coeff_ab_mo=np.random.rand(10, 10), olp_ab_ab=olp_ab_ab)
    _check_input(coeff_ab_mo=coeff_ab_mo, olp_ab_ab=olp_ab_ab)

    with pytest.raises(TypeError):
        _check_input(indices_span=indices_span.tolist())
    with pytest.raises(TypeError):
        _check_input(indices_span=indices_span.reshape(10, 1))
    with pytest.raises(TypeError):
        _check_input(indices_span=indices_span.astype(int))
    _check_input(indices_span=indices_span)

    with pytest.raises(ValueError):
        bad_indices_span = np.random.rand(9) < 0.5
        _check_input(indices_span=bad_indices_span, coeff_ab_mo=coeff_ab_mo)
    _check_input(indices_span=indices_span, coeff_ab_mo=coeff_ab_mo)


def test_make_mmo():
    """Test orbtools.quasi.make_mmo."""
    # olp_aao_ab
    olp_aao_ab = np.random.rand(5, 10)
    # coeff_ab_mo
    olp_ab_ab, *_ = np.linalg.svd(np.random.rand(10, 10))
    olp_ab_ab = (olp_ab_ab * np.random.rand(10)).dot(olp_ab_ab.T)
    norm_ab = np.diag(olp_ab_ab) ** (-0.5)
    olp_ab_ab *= norm_ab[:, None] * norm_ab[None, :]
    coeff_ab_mo = olp_ab_ab.dot(np.random.rand(10, 10))
    coeff_ab_mo *= np.diag(coeff_ab_mo.T.dot(olp_ab_ab).dot(coeff_ab_mo)) ** (-0.5)
    # indices_span
    indices_span = np.array([True] * 5 + [False] * 5)

    with pytest.raises(TypeError):
        make_mmo(olp_aao_ab, coeff_ab_mo, indices_span, dim_mmo=8.0)
    with pytest.raises(ValueError):
        make_mmo(olp_aao_ab, coeff_ab_mo, indices_span, dim_mmo=11)
    with pytest.raises(ValueError):
        make_mmo(olp_aao_ab, coeff_ab_mo, indices_span, dim_mmo=4)

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
    """Test orbtools.quasi.make_mmo against output from old code.

    Old code can be found in https://github.com/QuantumElephant/dumbo/tree/master/quasibasis.

    """
    # compare against reference generated using old code
    current_dir = os.path.dirname(__file__)
    olp_aao_ab = np.load(os.path.join(current_dir, "naclo4_olp_aao_ab.npy"))
    coeff_ab_mo = np.load(os.path.join(current_dir, "naclo4_coeff_ab_mo.npy"))
    occupations = np.load(os.path.join(current_dir, "naclo4_occupations.npy"))
    indices_span = occupations > 0

    coeff_ab_mmo = np.load(os.path.join(current_dir, "naclo4_coeff_ab_mmo.npy"))
    assert np.allclose(coeff_ab_mmo, make_mmo(olp_aao_ab, coeff_ab_mo, indices_span))


def test_quambo_old_code():
    """Test orbtools.quasi.quambo against output from old code.

    Old code can be found in https://github.com/QuantumElephant/dumbo/tree/master/quasibasis.

    """
    current_dir = os.path.dirname(__file__)
    coeff_ab_mo = np.load(os.path.join(current_dir, "naclo4_coeff_ab_mo.npy"))
    olp_ab_ab = np.load(os.path.join(current_dir, "naclo4_olp_ab_ab.npy"))
    olp_aao_ab = np.load(os.path.join(current_dir, "naclo4_olp_aao_ab.npy"))
    occupations = np.load(os.path.join(current_dir, "naclo4_occupations.npy"))
    indices_span = occupations > 0

    coeff_ab_quambo = np.load(os.path.join(current_dir, "naclo4_coeff_ab_quambo.npy"))
    olp_aao_ab = np.load(os.path.join(current_dir, "naclo4_olp_aao_ab.npy"))
    assert np.allclose(coeff_ab_quambo, quambo(olp_ab_ab, olp_aao_ab, coeff_ab_mo, indices_span))


def test_quao_old_code():
    """Test orbtools.quasi.quao against output from old code.

    Old code can be found in https://github.com/QuantumElephant/dumbo/tree/master/quasibasis.

    """
    current_dir = os.path.dirname(__file__)
    coeff_ab_mo = np.load(os.path.join(current_dir, "naclo4_coeff_ab_mo.npy"))
    olp_ab_ab = np.load(os.path.join(current_dir, "naclo4_olp_ab_ab.npy"))
    olp_aao_ab = np.load(os.path.join(current_dir, "naclo4_olp_aao_ab.npy"))
    olp_aao_aao = np.load(os.path.join(current_dir, "naclo4_olp_aao_aao.npy"))
    occupations = np.load(os.path.join(current_dir, "naclo4_occupations.npy"))
    indices_span = occupations > 0

    coeff_ab_quao = np.load(os.path.join(current_dir, "naclo4_coeff_ab_quao.npy"))
    olp_aao_ab = np.load(os.path.join(current_dir, "naclo4_olp_aao_ab.npy"))
    assert np.allclose(
        coeff_ab_quao, quao(olp_ab_ab, olp_aao_ab, olp_aao_aao, coeff_ab_mo, indices_span)
    )


def test_quambo():
    """Test orbtools.quasi.quambo against literature value for Mulliken populations.

    References
    ----------
    .. [1] Janowski, T. Near equivalence of intrinsic atomic orbitals and quasiatomic orbitals.
        JCTC, 2014, 10, 3085-3091.

    """
    current_dir = os.path.dirname(__file__)
    coeff_ab_mo = np.load(os.path.join(current_dir, "naclo4_coeff_ab_mo.npy"))
    olp_ab_ab = np.load(os.path.join(current_dir, "naclo4_olp_ab_ab.npy"))
    olp_aao_ab = np.load(os.path.join(current_dir, "naclo4_olp_aao_ab.npy"))
    occupations = np.load(os.path.join(current_dir, "naclo4_occupations.npy"))
    indices_span = occupations > 0
    ab_atom_indices = np.load(os.path.join(current_dir, "naclo4_qab_atom_indices.npy"))

    olp_ab_omo = olp_ab_ab.dot(coeff_ab_mo[:, indices_span])

    coeff_ab_quambo = quambo(olp_ab_ab, olp_aao_ab, coeff_ab_mo, indices_span)
    olp_quambo_quambo = coeff_ab_quambo.T.dot(olp_ab_ab).dot(coeff_ab_quambo)
    olp_quambo_omo = coeff_ab_quambo.T.dot(olp_ab_omo)
    coeff_quambo_omo = project(olp_quambo_quambo, olp_quambo_omo)
    pop = mulliken_populations(
        coeff_quambo_omo, occupations[indices_span], olp_quambo_quambo, 6, ab_atom_indices
    )
    partial_pop = np.array([11, 17, 8, 8, 8, 8]) - pop
    assert np.allclose(
        partial_pop, np.array([0.968, 2.454, -0.807, -0.904, -0.904, -0.807]), atol=1e-3
    )


def test_quao():
    """Test orbtools.quasi.quambo against literature value for Mulliken populations.


    References
    ----------
    .. [1] Janowski, T. Near equivalence of intrinsic atomic orbitals and quasiatomic orbitals.
        JCTC, 2014, 10, 3085-3091.

    """
    current_dir = os.path.dirname(__file__)
    coeff_ab_mo = np.load(os.path.join(current_dir, "naclo4_coeff_ab_mo.npy"))
    olp_ab_ab = np.load(os.path.join(current_dir, "naclo4_olp_ab_ab.npy"))
    olp_aao_ab = np.load(os.path.join(current_dir, "naclo4_olp_aao_ab.npy"))
    olp_aao_aao = np.load(os.path.join(current_dir, "naclo4_olp_aao_aao.npy"))
    occupations = np.load(os.path.join(current_dir, "naclo4_occupations.npy"))
    indices_span = occupations > 0
    ab_atom_indices = np.load(os.path.join(current_dir, "naclo4_qab_atom_indices.npy"))

    olp_ab_omo = olp_ab_ab.dot(coeff_ab_mo[:, indices_span])

    coeff_ab_quao = quao(olp_ab_ab, olp_aao_ab, olp_aao_aao, coeff_ab_mo, indices_span)
    olp_quao_quao = coeff_ab_quao.T.dot(olp_ab_ab).dot(coeff_ab_quao)
    olp_quao_omo = coeff_ab_quao.T.dot(olp_ab_omo)
    coeff_quao_omo = project(olp_quao_quao, olp_quao_omo)
    pop = mulliken_populations(
        coeff_quao_omo, occupations[indices_span], olp_quao_quao, 6, ab_atom_indices
    )
    partial_pop = np.array([11, 17, 8, 8, 8, 8]) - pop
    assert np.allclose(
        partial_pop, np.array([0.967, 2.498, -0.819, -0.914, -0.914, -0.819]), atol=1e-3
    )
