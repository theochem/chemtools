"""Tests for orbtools.orthogonalization."""
import chemtools.orbtools.orthogonalization as orth
import numpy as np
from numpy.testing import assert_raises


def test_eigh():
    """Test orbtools.orthogonalization.eigh."""
    assert_raises(TypeError, orth.eigh, np.random.rand(3, 3).tolist())
    assert_raises(TypeError, orth.eigh, np.random.rand(3, 3, 3))
    assert_raises(ValueError, orth.eigh, np.random.rand(3, 5))
    assert_raises(ValueError, orth.eigh, np.random.rand(3, 3))
    matrix = np.random.rand(3, 3)
    assert_raises(TypeError, orth.eigh, matrix + matrix.T, threshold=None)
    matrix = np.random.rand(3, 3)
    assert_raises(ValueError, orth.eigh, matrix + matrix.T, threshold=-2)
    # positive semidefinite
    matrix = np.arange(9).reshape(3, 3)
    matrix = matrix.dot(matrix.T)
    eigval, eigvec = orth.eigh(matrix)
    assert np.allclose(eigval, np.array([2.02399203e02, 1.60079682e00]))
    assert np.allclose(
        eigvec,
        np.array(
            [[-0.13511895, -0.90281571], [-0.49633514, -0.29493179], [-0.85755134, 0.31295213]]
        ),
    )
    # reconstruct random symmetric positive semidefinite
    matrix = np.random.rand(100, 100)
    matrix = matrix.dot(matrix.T)
    eigval, eigvec = orth.eigh(matrix)
    assert np.allclose((eigvec * eigval).dot(eigvec.T), matrix)
    # reconstruct random symmetric
    matrix = np.random.rand(100, 100)
    matrix = matrix + matrix.T
    eigval, eigvec = orth.eigh(matrix)
    assert np.allclose((eigvec * eigval).dot(eigvec.T), matrix)


def test_svd():
    """Test orbtools.orthogonalization.svd."""
    assert_raises(TypeError, orth.svd, np.random.rand(3, 3).tolist())
    assert_raises(TypeError, orth.svd, np.random.rand(3, 3, 3))
    matrix = np.random.rand(10, 10)
    u, sigma, vdagger = orth.svd(matrix)
    assert np.allclose((u * sigma).dot(vdagger), matrix)

    matrix = np.array(
        [
            [0.20090975, 0.97054546, 0.28661335, 0.07573257, 0.3917035, 0.75842177],
            [0.87478886, 0.95606901, 0.39057639, 0.36235354, 0.1355327, 0.15858099],
            [0.02462707, 0.05013185, 0.00880004, 0.13152568, 0.9655135, 0.79478605],
            [0.97382223, 0.27750593, 0.64722406, 0.60750097, 0.78990422, 0.26197507],
        ]
    )
    u, sigma, vdagger = orth.svd(matrix, threshold=1)
    assert u.shape == (4, 2)
    assert np.allclose(sigma, np.array([2.36256353, 1.17763142]))
    assert vdagger.shape == (2, 6)


def test_power_symmetric():
    """Test orbtools.orthogonalization.power_symmetric."""
    matrix = np.random.rand(5, 5)
    # positive semidefinite
    matrix = matrix.dot(matrix.T)
    assert np.allclose(orth.power_symmetric(matrix, 2), matrix.dot(matrix))
    assert np.allclose(orth.power_symmetric(matrix, -1).dot(matrix), np.identity(5))
    power_matrix = orth.power_symmetric(matrix, 0.5)
    assert np.allclose(power_matrix.dot(power_matrix), matrix)
    power_matrix = orth.power_symmetric(matrix, 1.0 / 3.0)
    assert np.allclose(power_matrix.dot(power_matrix).dot(power_matrix), matrix)
    # not positive semidifinite (probably)
    matrix = matrix + matrix.T
    assert np.allclose(orth.power_symmetric(matrix, 2), matrix.dot(matrix))
    assert np.allclose(orth.power_symmetric(matrix, -1).dot(matrix), np.identity(5))
    # following crashes because it has negative eigenvalues
    matrix = np.random.rand(100, 100)
    matrix = matrix + matrix.T
    assert_raises(ValueError, orth.power_symmetric, matrix, 0.5)
