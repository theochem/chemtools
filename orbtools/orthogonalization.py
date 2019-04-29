"""This code has been copied from https://github.com/QuantumElephant/dumbo."""


class OrthogonalizationError(Exception):
    """ Exception class for errors in the orthogonalization module

    """

    pass


def eigh(matrix, threshold=1e-9):
    """Returns eigenvalues and eigenvectors of a Hermitian matrix where the
    eigenvectors (and eigenvalues) with eigenvalues less than the threshold are
    removed.

    Parameters
    ----------
    matrix : np.ndarray(N,N)
        Square Hermitian matrix
    threshold : {1e-9, float}
        Eigenvalues (and corresponding eigenvectors) below this threshold are discarded

    Returns
    -------
    eigval : np.ndarray(K,)
        Eigenvalues sorted in decreasing order
    eigvec : np.ndarray(N,K)
        Matrix where the columns are the corresponding eigenvectors to the eigval

    Raises
    ------
    OrthogonalizationError
        If matrix is not a square two dimensional numpy array

    NOTE
    ----
    This code mainly uses numpy.eigh
    """
    if not isinstance(matrix, np.ndarray):
        raise OrthogonalizationError("Unsupported matrix type, {0}".format(type(matrix)))
    if len(matrix.shape) != 2 or matrix.shape[0] != matrix.shape[1]:
        raise OrthogonalizationError("Unsupported matrix shape, {0}".format(matrix.shape))
    eigval, eigvec = np.linalg.eigh(matrix)
    # discard eigenvalues less than threshold
    kept_indices = np.abs(eigval) > threshold
    if np.sum(-kept_indices) > 0:
        print(
            "WARNING: Discarded {0} eigenvalues (threshold= {1}):\n"
            "{2}".format(sum(-kept_indices), threshold, eigval[-kept_indices])
        )
    if np.any(eigval < -threshold):
        print(
            "WARNING: {0} eigenvalues are quite negative:\n"
            "{1}".format(np.sum(eigval < -threshold), eigval[eigval < -threshold])
        )
    eigval, eigvec = eigval[kept_indices], eigvec[:, kept_indices]
    # sort it by decreasing eigenvalue
    sorted_indices = np.argsort(eigval, kind="quicksort")[::-1]
    return eigval[sorted_indices], eigvec[:, sorted_indices]


def svd(matrix, threshold=1e-9):
    """Returns singular values and vectors with singular values greater than
    threshold

    Parameters
    ----------
    matrix : np.ndarray(N,M)
        Matrix (not necessarily square)

    threshold : {1e-9, float}
        Singularvalues (and corresponding singularvectors) below this threshold are
        discarded

    Returns
    -------
    u : np.ndarray(N,K)
        Left singular matrix that is unitary
    sigma : np.ndarray(K,)
        Singular values sorted in decreasing order
    vdagger : np.ndarray(K,M)
        Right singular matrix that is unitary

    Raises
    ------
    OrthogonalizationError
        If matrix is not a two dimensional numpy array

    NOTE
    ----
    This code uses numpy.linalg.svd

    """
    if not isinstance(matrix, np.ndarray):
        raise OrthogonalizationError("Unsupported matrix type, {0}".format(type(matrix)))
    if len(matrix.shape) != 2:
        raise OrthogonalizationError("Unsupported matrix shape, {0}".format(matrix.shape))
    u, sigma, vdagger = np.linalg.svd(matrix, full_matrices=False)
    # discard eigenvalues less than threshold
    kept_indices = sigma > threshold
    if np.sum(-kept_indices) > 0:
        print(
            "WARNING: Discarded {0} eigenvalues (threshold= {1}):\n"
            "{2}".format(sum(-kept_indices), threshold, sigma[-kept_indices])
        )
    u, sigma, vdagger = u[:, kept_indices], sigma[kept_indices], vdagger[kept_indices, :]
    # sort it by decreasing singular
    sorted_indices = np.argsort(sigma, kind="quicksort")[::-1]
    return u[:, sorted_indices], sigma[sorted_indices], vdagger[sorted_indices, :]


def power_symmetric(matrix, k, threshold_eig=1e-9, threshold_symm=1e-10):
    """ Return kth power of a symmetric matrix

    Parameters
    ----------
    matrix : np.ndarray(N,N)
        Symmetric matrix
    k : float
        The power of the matrix
    threshold_eig : {1e-9, float}
        In the eigenvalue decomposition, the eigenvalues (and corresponding
        eigenvectors) that are less than the threshold are discarded
    threshold_symm : {1e-10, float}
        Used to check that the matrix is symmetric.

    Returns
    -------
    answer : np.ndarray(N,N)
        The matrix raised to the kth power

    Raises
    ------
    OrthogonalizationError
        If the matrix is not a square two dimensional numpy array
        If the maximum of the absolute difference of the matrix with its transpose
        is greater than the threshold_symm, an error is raised
        If the power is a fraction and the any eigenvalues are negative
    """
    if not isinstance(matrix, np.ndarray):
        raise OrthogonalizationError("Unsupported matrix type, {0}".format(type(matrix)))
    if len(matrix.shape) != 2 or matrix.shape[0] != matrix.shape[1]:
        raise OrthogonalizationError("Unsupported matrix shape, {0}".format(matrix.shape))
    if np.max(np.abs(matrix.T - matrix)) > threshold_symm:
        raise OrthogonalizationError("Matrix is not symmetric")
    eigval, eigvec = eigh(matrix, threshold=threshold_eig)
    if k % 1 != 0 and np.any(eigval < 0):
        raise OrthogonalizationError(
            "Fractional power of negative eigenvalues. " "Imaginary numbers not supported."
        )
    return (eigvec * (eigval ** k)).dot(eigvec.T)
