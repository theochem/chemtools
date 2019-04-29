"""Mulliken population analysis.

This code has been copied from https://github.com/QuantumElephant/dumbo.

"""
import numpy as np


class Mulliken:
    r"""Class for applying Mulliken analysis.

    ..math::
        \mo_i = \sum_i \ketf{\ab_j} C_{ji}

    ..math::
        1 &= \braket{\mo_i}{\mo_i}\\
        N &= \sum_i^{occ} n_i \braket{\mo_i}{\mo_i}\\

    where N is the number of electrons and :math:`n_i` is the occupation number of
    molecular orbital i.

    ..math::
        N &= \sum_i^{occ} n_i \braket{\mo_i}{\mo_i}\\
        &= \sum_i^{occ} n_i \sum_{jk} C_{ij}^\dagger \braket{\ab_j|\ab_k} C_{ki}\\
        &= \sum_{jk} \braket{\ab_j}{\ab_k} \sum_i^{occ} C_{ki} n_i C_{ij}^\dagger\\
        &= \sum_{jk} S_{jk} X_{kj}\\

    where

    ..math::
        S_{jk} &= \braket{\ab_k}{\ab_j}\\
        X_{kj} &= \sum_i^{occ} C_{ki} n_i C_{ij}^\dagger \\

    Divide up the sum of N into corresponding contribution from each atom

    ..math::
        N = \sum_A \sum_{j,k} a_{jk}^A S_{jk} X_{kj}

    where (j,k) and :math:`a_{jk}^A` corresponds to atom A and depends on the
    ``partition'' scheme of :math:`S_{\ab}`.

    Note that
    ..math::
        N &= \sum_A \sum_{jk} a_{jk}^A S_{jk} X_{kj}\\
        &= \sum_{jk} (\sum_A a_{jk}^A) S_{jk} X_{kj}\\
        &= \sum_{jk} S_{jk} X_{kj}

    So :math:`\sum_A a_{jk}^A = 1` for any j and k.

    """
    def __init__(self, coeff_ab_mo, occupations, olp_ab_ab, num_center, basis_map, weights=None):
        """Initialize Mulliken analysis.

        Parameters
        ----------
        coeff_ab_mo : np.ndarray(K, N)
            Transformation matrix from the atomic basis to molecular orbitals
        occupations : np.ndarray(N, )
            Occupation numbers of each of the molecular orbitals
        olp_ab_ab : np.ndarray(K, K)
            Overlap matrix of the molecular orbitals
        weights : {None, np.ndarray(N, N)}

        Raises
        ------
        TypeError
            If `coeff_ab_mo` is not a numpy array.
            If `coeff_ab_mo` is not a two-dimensional numpy array.
            If `occupations` is not a numpy array.
            If `occupations` is not a one-dimensional numpy array.
            If `olp_ab_ab` is not a numpy array.
            If `olp_ab_ab` is not a two-dimensional numpy array.
            If number of atom centers is not an integer
            If the basis mapping is not iterable
            If `weights` is not iterable.
            If `weights` contains objects that are not numpy arrays.
            If `weights` is not a two dimensional array
        ValueError
            If `olp_ab_ab` is not symmetric.
            If `olp_ab_ab` is not normalized.
            If molecular orbitals are not normalized
            If overlap matrix is not square.
            If number of atomic orbitals in the transformation matrix and overlap matrix are not
            equal.
            If number of molecular orbitals in the transformation matrix and occupations are not
            equal.
            If there are negative occupation numbers
            If the indices used in the basis mapping must be within the number of atom centers
            If `weights` does not have an entry for each atom center
            If `weights` does not have the same shape as the overlap matrix
            If `weights` is not symmeteric
            If `weights` is not normalized (add up to 1)
        """
        # Check if numpy arrays of right shape
        if not isinstance(coeff_ab_mo, np.ndarray):
            raise TypeError('`coeff_ab_mo` is not a numpy array.')
        elif not len(coeff_ab_mo.shape) == 2:
            raise TypeError('`coeff_ab_mo` is not a two-dimensional numpy array.')

        if not isinstance(occupations, np.ndarray):
            raise TypeError('`occupations` is not a numpy array.')
        elif not len(occupations.shape) == 1:
            raise TypeError('`occupations` is not a one-dimensional numpy array.')

        if not isinstance(olp_ab_ab, np.ndarray):
            raise TypeError('`olp_ab_ab` is not a numpy array.')
        elif not len(olp_ab_ab.shape) == 2:
            raise TypeError('`olp_ab_ab` is not a two-dimensional numpy array.')

        # Check if overlap is symmetric
        if not np.allclose(olp_ab_ab, olp_ab_ab.T):
            raise ValueError('`olp_ab_ab` is not symmetric.')
        # Check if overlap is normalized
        if not np.allclose(np.diag(olp_ab_ab), np.ones(coeff_ab_mo.shape[0])):
            raise ValueError('`olp_ab_ab` is not normalized.')
        # Check if molecular orbitals are normalized
        if not np.allclose(np.diag(coeff_ab_mo.T.dot(olp_ab_ab).dot(coeff_ab_mo)), 1):
            raise ValueError('Molecular orbitals are not normalized.')
        # Check if given values are consistent
        if not olp_ab_ab.shape[0] == olp_ab_ab.shape[1]:
            raise ValueError('Overlap matrix is not square.')
        if not coeff_ab_mo.shape[0] == olp_ab_ab.shape[0]:
            raise ValueError('Number of atomic orbitals in the transformation matrix and overlap '
                             'matrix are not equal.')
        if not coeff_ab_mo.shape[1] == occupations.size:
            raise ValueError('Number of molecular orbitals in the transformation matrix and '
                             'occupations are not equal.')
        # Check if occupations look okay
        if not np.all(occupations >= 0):
            raise ValueError('Negative occupation numbers.')
        if np.any(occupations > 2):
            print('Warning: Atleast one occupation number exceeds 2.')
        # Check number of center
        if not isinstance(num_center, int):
            raise TypeError('Number of atom centers must be an integer')
        # Check basis mapping
        if not hasattr(basis_map, '__iter__'):
            raise TypeError('Basis mapping must be iterable.')
        if not all(i in range(num_center) for i in basis_map):
            raise ValueError('Indices used in the basis mapping must be within the number of atom'
                             'centers')

        if weights is not None:
            if not hasattr(weights, '__iter__'):
                raise TypeError('`weights` must be iterable.')
            if not len(weights) == num_center:
                raise ValueError('`weights` must have one entry for each atom center.')
            for weight in weights:
                if not isinstance(weight, np.ndarray):
                    raise TypeError('`weights` must contain numpy arrays.')
                if not len(weight.shape) == 2:
                    raise TypeError('`weights` must contain two-dimensional numpy arrays.')
                if weight.shape != olp_ab_ab.shape:
                    raise ValueError('`weights` must have the same shape as the overlap matrix.')
                if not np.allclose(weight, weight.T):
                    raise ValueError('`weights` must be symmetric.')
            # check if normalized
            if not np.allclose(sum(weights), np.ones(olp_ab_ab.shape)):
                raise ValueError('`weights` must be normalized.')

        # Initialize
        #  Remove unoccupied
        ind_occ = occupations > 0
        self.coeff_ab_mo = coeff_ab_mo[:, ind_occ]
        self.occupations = occupations[ind_occ]
        self.olp_ab_ab = olp_ab_ab
        self.num_center = num_center
        self.basis_map = tuple(basis_map)

        if weights is None:
            self.weights = self.weights_default(num_center, basis_map)
        else:
            self.weights = weights

    @property
    def num_ab(self):
        return self.coeff_ab_mo.shape[0]

    def weights_default(self, num_center, basis_map):
        r"""Return the default weights for Mulliken analysis.

        ..math::
            a_{jk}^A &= 1 \mbox{ if $j \in A$ and $k \in A$}\\
            &= 0.5 \mbox{ if $j \in A$ and $k \not\in A$}\\
            &= 0.5 \mbox{ if $j \not\in A$ and $k \in A$}\\
            &= 0 \mbox{ if $j \not\in A$ and $k \not\in A$\\

        Parameters
        ----------
        num_center : int
            Number of atom centers.
        basis_map : iterable
            List of indices of the center of which the corresponding basis function (atomic basis)
            belong

        Returns
        -------
        weights : list of np.ndarray
            Population weight matrix for each atom center

        Raises
        ------
        TypeError
            If `basis_map` is not iterable
            If there are indices in `basis_map` that are outside the range of `num_center`

        Example
        -------
        If the 1st through 4th basis functions belongs on the 1st atom, and
        the 5th through 10th basis functions belong on the 2nd atom, then
        basis_map = [0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
        """
        output = []
        if not hasattr(basis_map, '__iter__'):
            raise TypeError('`basis_map` must be iterable')
        if not all(i in range(num_center) for i in basis_map):
            raise TypeError('Each index in `basis_map` must be within range of `num_center`')
        for ind_center in range(num_center):
            weight = np.zeros([self.num_ab]*2)
            indices = np.array([i for i, j in enumerate(basis_map) if j == ind_center])
            if indices.size != 0:
                weight[indices, :] += 0.5
                weight[:, indices] += 0.5
            output.append(weight)
        return output

    @staticmethod
    def weights_olp(num_center, aao_map, olp_aao_basis):
        """Return the weights modified using the overlap matrix.

        Weights are distributed by the overlap of the given orbital with the AAO's. If the overlap
        is large, then a larger weight is attributed.

        Parameters
        ----------
        num_center : int
            Number of atom centers.
        aao_map : iterable
            List of indices of the center of which the corresponding basis function (AAO) belong
        olp_aao_basis : np.ndarray
            Overlap matrix between the AAO and the AB

        Raises
        ------
        TypeError
            If `aao_map` is not an iterable
        ValueError
            If `aao_map` does not have indices that are within range of `num_center`
        """
        if not hasattr(aao_map, '__iter__'):
            raise TypeError('`aao_map` is not an iterable')
        if not all(i in range(num_center) for i in aao_map):
            raise ValueError('`aao_map` does not have indices that are within range of '
                             '`num_center`')
        output = []
        num_basis = olp_aao_basis.shape[1]
        olp_aao_basis = np.abs(olp_aao_basis)**2
        total = np.sum(olp_aao_basis, axis=0)
        total = total.reshape(1, total.size)
        total_2d = total+total.T
        total_2d[np.abs(total_2d) < 1e-9] = 1
        for ind_center in range(num_center):
            indices = np.array([i for i, j in enumerate(aao_map) if j == ind_center])
            weight = np.zeros([num_basis]*2)
            if indices.size != 0:
                top = np.sum(olp_aao_basis[indices, :], axis=0)
                top = top.reshape(1, top.size)
                top_2d = top+top.T
                weight += top_2d/total_2d
            output.append(weight)
        hack = sum(output)
        for i in range(num_center):
            output[i][hack < 1e-9] = 1.0/num_center
        return output

    @staticmethod
    def weights_chi(num_center, basis_map, electronegs):
        """Return the weights modified using the electronegativity.

        Weights are distributed by the electronegativity of the given atom center. If the atom is
        electronegative, then a larger weight is attributed.

        Parameters
        ----------
        num_center : int
            Number of atom centers.
        basis_map : iterable
            List of indices of the center of which the corresponding basis function (AAO) belong
        electronegs : np.ndarray
            Electronegativity of the atom centers

        Raises
        ------
        TypeError
            If `aao_map` is not an iterable
        ValueError
            If `aao_map` does not have indices that are within range of `num_center`
        """
        if not hasattr(basis_map, '__iter__'):
            raise TypeError('`basis_map` is not an iterable')
        if not all(i in range(num_center) for i in basis_map):
            raise ValueError('`basis_map` does not have indices that are within range of '
                             '`num_center`')
        output = []
        num_basis = len(basis_map)
        electronegs = electronegs.reshape(1, electronegs.size)
        total_electronegs = electronegs+electronegs.T
        for ind_center in range(num_center):
            weight = np.zeros([num_basis]*2)
            indices = np.array([i for i, j in enumerate(basis_map) if j == ind_center])
            if indices.size != 0:
                weight[indices, :] += electronegs.T[indices, :]/total_electronegs[indices, :]
                weight[:, indices] += electronegs[:, indices]/total_electronegs[:, indices]
            output.append(weight)
        return output

    @staticmethod
    def weights_olp_chi(num_center, aao_map, olp_aao_basis, electronegs):
        """Return the weights modified using the overlap matrix and the electronegativities.

        Weights are distributed by the overlap of the given orbital with the AAO's and the
        electronegativities of the atoms. If the overlap and the electronegativity is large, then a
        larger weight is attributed.

        Parameters
        ----------
        num_center : int
            Number of atom centers.
        aao_map : iterable
            List of indices of the center of which the corresponding basis function (AAO) belong
        olp_aao_basis : np.ndarray
            Overlap matrix between the AAO and the AB
        electronegs : np.ndarray
            Electronegativity of the atom centers

        Raises
        ------
        TypeError
            If `aao_map` is not an iterable
        ValueError
            If `aao_map` does not have indices that are within range of `num_center`
        """
        if not hasattr(aao_map, '__iter__'):
            raise TypeError('`aao_map` is not an iterable')
        if not all(i in range(num_center) for i in aao_map):
            raise ValueError('`aao_map` does not have indices that are within range of '
                             '`num_center`')
        output = []
        num_basis = olp_aao_basis.shape[1]
        electronegs = electronegs.reshape(electronegs.size, 1)
        olp_aao_basis = olp_aao_basis**2*electronegs
        total = np.sum(olp_aao_basis, axis=0)
        total = total.reshape(1, total.size)
        total_2d = total+total.T
        total_2d[np.abs(total_2d) < 1e-9] = 1
        for ind_center in range(num_center):
            indices = np.array([i for i, j in enumerate(aao_map) if j == ind_center])
            weight = np.zeros([num_basis]*2)
            if indices.size != 0:
                top = np.sum(olp_aao_basis[indices, :], axis=0)
                top = top.reshape(1, top.size)
                top_2d = top+top.T
                weight += top_2d/total_2d
            output.append(weight)
        hack = sum(output)
        for i in range(num_center):
            output[i][hack < 1e-9] = 1.0/num_center
        return output

    def get_population(self):
        """ Returns the numpy array for the atomic populations.

        Returns
        -------
        populations : np.ndarray
            One dimension numpy array that contains the mulliken charges or partial charges of each
            center

        Raises
        ------
        ValueError
            If the total population is not equal to the total number of electrons
        """
        output = []
        num_elec = np.sum(self.occupations)
        X = (self.coeff_ab_mo*self.occupations).dot(self.coeff_ab_mo.T)
        for weight in self.weights:
            output.append(np.sum(self.olp_ab_ab*X.T*weight))
        if not abs(num_elec-sum(output)) < 1e-6:
            raise ValueError('Population does not match up with the number of electrons.')
        return np.array(output)
