"""Tools for matrix decomposition and power."""
import numpy as np


def eigh(matrix, threshold=1e-9):
    """Return the eigenvalues and eigenvectors of a Hermitian matrix.

    Eigenvalues whose absolute values are less than the threshold are discarded as well as the
    corresponding eigenvectors.

    Parameters
    ----------
    matrix : np.ndarray(N, N)
        Square Hermitian matrix.
    threshold : {1e-9, float}
        Eigenvalues (and corresponding eigenvectors) below this threshold are discarded.

    Returns
    -------
    eigval : np.ndarray(K,)
        Eigenvalues sorted in decreasing order.
    eigvec : np.ndarray(N,K)
        Matrix where the columns are the corresponding eigenvectors to the eigval.

    Raises
    ------
    TypeError
        If `matrix` is not a two-dimensional numpy array.
        If `threshold` is not an integer or a float.
    ValueError
        If `matrix` is not a square matrix.
        If `matrix` is not Hermitian.
        If `threshold` is negative.

    Warns
    -----
    If there are any negative eigenvalues (beyond the given threshold).
    If any eigenvalues and eigenvectors are discarded.

    Note
    ----
    This code mainly uses numpy.eigh

    """
    if not (isinstance(matrix, np.ndarray) and matrix.ndim == 2):
        raise TypeError("Given matrix must be a two-dimensional numpy array.")
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Given matrix must be square.")
    if not np.allclose(matrix.conjugate().T, matrix):
        raise ValueError("Given matrix must be Hermitian.")
    if not isinstance(threshold, (int, float)):
        raise TypeError("Given threshold must be an integer or a float.")
    if threshold < 0:
        raise ValueError("Given threshold must be positive.")

    eigval, eigvec = np.linalg.eigh(matrix)
    # NOTE: it is assumed that the np.linalg.eigh sorts the eigenvalues in increasing order.

    neg_indices = eigval < -threshold
    if np.any(neg_indices):
        print(
            "WARNING: {0} eigenvalues are negative (less than the threshold {1}):\n"
            "{2}".format(np.sum(neg_indices), -threshold, eigval[neg_indices])
        )

    kept_indices = np.abs(eigval) > threshold
    if np.sum(~kept_indices) > 0:
        print(
            "WARNING: Discarded {0} eigenvalues because they are less than the threshold {1}:\n"
            "{2}".format(np.sum(~kept_indices), threshold, eigval[~kept_indices])
        )

    eigval, eigvec = eigval[kept_indices], eigvec[:, kept_indices]
    return eigval[::-1], eigvec[:, ::-1]


def svd(matrix, threshold=1e-9):
    """Return the singular values and singular vectors of the given matrix.

    Singular values whose absolute values are less than the threshold are discarded as well as the
    corresponding singular vectors.


    Parameters
    ----------
    matrix : np.ndarray(N, M)
        Matrix.
    threshold : {1e-9, float}
        Singular values (and corresponding singular vectors) below this threshold are discarded.

    Returns
    -------
    u : np.ndarray(N, K)
        Left singular matrix.
    sigma : np.ndarray(K,)
        Singular values sorted in decreasing order.
    vdagger : np.ndarray(K, M)
        Right singular matrix.

    Raises
    ------
    TypeError
        If `matrix` is not a two-dimensional numpy array.

    Warns
    -----
    If any singular values and singular vectors are discarded.

    Note
    ----
    This code uses numpy.linalg.svd

    """
    # pylint: disable=C0103
    if not (isinstance(matrix, np.ndarray) and matrix.ndim == 2):
        raise TypeError("Given matrix must be a two-dimensional numpy array.")

    u, sigma, vdagger = np.linalg.svd(matrix, full_matrices=False)
    # NOTE: it is assumed that the np.linalg.svd sorts the singular values in descending order.

    kept_indices = sigma > threshold
    if np.sum(~kept_indices) > 0:
        print(
            "WARNING: Discarded {0} singular values because they are less than the threshold {1}:\n"
            "{2}".format(np.sum(~kept_indices), threshold, sigma[~kept_indices])
        )

    u, sigma, vdagger = u[:, kept_indices], sigma[kept_indices], vdagger[kept_indices, :]

    return u, sigma, vdagger


def power_symmetric(matrix, k, threshold=1e-9):
    """Return the kth power of the given symmetric matrix.

    Parameters
    ----------
    matrix : np.ndarray(N, N)
        Symmetric matrix.
    k : {int, float}
        Power of the matrix.
    threshold : {1e-9, float}
        In the eigenvalue decomposition, the eigenvalues (and corresponding eigenvectors) that are
        less than the threshold are discarded.

    Returns
    -------
    matrix_power : np.ndarray(N, N)
        Matrix raised to the kth power.

    Raises
    ------
    ValueError
        If the `k` is a fraction and matrix has negative eigenvalues.

    """
    eigval, eigvec = eigh(matrix, threshold=threshold)
    if k % 1 != 0 and np.any(eigval < 0):
        raise ValueError(
            "Given matrix has negative eigenvalues. Fractional powers of negative eigenvalues are "
            "not supported."
        )
    return (eigvec * (eigval ** k)).dot(eigvec.T)
