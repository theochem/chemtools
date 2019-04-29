import numpy as np
import quasibasis.orthogonalization as orth


def test_eigh():
    # positive semidefinite
    A = np.arange(9).reshape(3, 3)
    B = A.dot(A.T)
    eigval, eigvec = orth.eigh(B)
    assert np.allclose(eigval, np.array([2.02399203e02, 1.60079682e00]))
    assert np.allclose(
        eigvec,
        np.array(
            [[-0.13511895, -0.90281571], [-0.49633514, -0.29493179], [-0.85755134, 0.31295213]]
        ),
    )
    # reconstruct random symmetric positive semidefinite
    A = np.random.rand(100, 100)
    B = A.dot(A.T)
    eigval, eigvec = orth.eigh(B)
    assert np.allclose((eigvec * eigval).dot(eigvec.T), B)
    # not positive semidefinite
    A = np.arange(9).reshape(3, 3)
    B = (A + A.T) / 2
    eigval, eigvec = orth.eigh(B)
    # reconstruct random symmetric positive semidefinite
    A = np.random.rand(100, 100)
    B = (A + A.T) / 2
    eigval, eigvec = orth.eigh(B)
    assert np.allclose((eigvec * eigval).dot(eigvec.T), B)


def test_svd():
    A = np.array(
        [
            [0.20090975, 0.97054546, 0.28661335, 0.07573257, 0.3917035, 0.75842177],
            [0.87478886, 0.95606901, 0.39057639, 0.36235354, 0.1355327, 0.15858099],
            [0.02462707, 0.05013185, 0.00880004, 0.13152568, 0.9655135, 0.79478605],
            [0.97382223, 0.27750593, 0.64722406, 0.60750097, 0.78990422, 0.26197507],
        ]
    )
    U, S, Vdag = orth.svd(A, threshold=1)
    assert np.allclose(
        U,
        np.array(
            [
                [-0.48550089, -0.16348557],
                [-0.52481293, 0.5508918],
                [-0.34346377, -0.81060264],
                [-0.60900978, 0.11275662],
            ]
        ),
    )
    assert np.allclose(S, np.array([2.36256353, 1.17763142]))
    assert np.allclose(
        Vdag,
        np.array(
            [
                [-0.49021672, 0.45762222],
                [-0.49064516, 0.30457239],
                [-0.31377732, 0.1988344],
                [-0.27177445, 0.12662799],
                [-0.45458249, -0.57993941],
                [-0.37415518, -0.55309861],
            ]
        ).T,
    )


def test_power_symmetric():
    A = np.random.rand(5, 5)
    # positive semidefinite
    B = A.dot(A.T)
    assert np.allclose(orth.power_symmetric(B, 2), B.dot(B))
    assert np.allclose(orth.power_symmetric(B, -1).dot(B), np.identity(5))
    C = orth.power_symmetric(B, 0.5)
    assert np.allclose(C.dot(C), B)
    C = orth.power_symmetric(B, 1.0 / 3.0)
    assert np.allclose(C.dot(C).dot(C), B)
    # not positive semidifinite (probably)
    B = (A + A.T) / 2
    assert np.allclose(orth.power_symmetric(B, 2), B.dot(B))
    assert np.allclose(orth.power_symmetric(B, -1).dot(B), np.identity(5))
    # following crashes because it has negative eigenvalues
    # C = orth.power_symmetric(B, 0.5)
    # assert np.allclose(C.dot(C), B)


test_eigh()
test_svd()
test_power_symmetric()
