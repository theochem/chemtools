"""Test orbtools.mulliken."""
import os

import numpy as np
from orbtools.mulliken import mulliken_populations, mulliken_populations_newbasis
from orbtools.quasi import project
import pytest


def test_mulliken_populations_input():
    """Test input checks in the orbtools.mulliken.mulliken_populations."""
    # get random unitary matrix
    unitary = np.linalg.svd(np.random.rand(20, 20))[0]
    # get random olp_ab_ab
    olp_ab_ab = (unitary * np.random.rand(20)).dot(unitary.T)
    norm = np.diag(olp_ab_ab) ** (-0.5)
    olp_ab_ab *= norm[:, None]
    olp_ab_ab *= norm[None, :]
    # get random mo's
    coeff_ab_mo = np.random.rand(20, 15) - 0.5
    coeff_ab_mo *= np.diag(coeff_ab_mo.T.dot(olp_ab_ab).dot(coeff_ab_mo)) ** (-0.5)
    occupations = np.random.rand(15)
    num_atoms = 4
    ab_atom_indices = np.array([0, 1, 2, 1, 1, 0, 2, 1, 0, 2, 1, 2, 0, 1, 2, 0, 3, 3, 1, 0])
    atom_weights = np.random.rand(4, 20, 20)
    atom_weights += np.swapaxes(atom_weights, 1, 2)
    atom_weights /= np.sum(atom_weights, axis=0)[None, :, :]

    with pytest.raises(TypeError):
        mulliken_populations(
            coeff_ab_mo.tolist(), occupations, olp_ab_ab, num_atoms, ab_atom_indices
        )
    with pytest.raises(TypeError):
        mulliken_populations(
            coeff_ab_mo.ravel(), occupations, olp_ab_ab, num_atoms, ab_atom_indices
        )
    with pytest.raises(TypeError):
        mulliken_populations(
            coeff_ab_mo.astype(int), occupations, olp_ab_ab, num_atoms, ab_atom_indices
        )

    with pytest.raises(TypeError):
        mulliken_populations(
            coeff_ab_mo, occupations.tolist(), olp_ab_ab, num_atoms, ab_atom_indices
        )
    with pytest.raises(TypeError):
        mulliken_populations(
            coeff_ab_mo, occupations[:, None], olp_ab_ab, num_atoms, ab_atom_indices
        )
    with pytest.raises(TypeError):
        mulliken_populations(
            coeff_ab_mo, occupations.astype(bool), olp_ab_ab, num_atoms, ab_atom_indices
        )

    with pytest.raises(TypeError):
        mulliken_populations(
            coeff_ab_mo, occupations, olp_ab_ab.tolist(), num_atoms, ab_atom_indices
        )
    with pytest.raises(TypeError):
        mulliken_populations(
            coeff_ab_mo, occupations, olp_ab_ab.ravel(), num_atoms, ab_atom_indices
        )
    with pytest.raises(TypeError):
        mulliken_populations(
            coeff_ab_mo, occupations, olp_ab_ab.astype(int), num_atoms, ab_atom_indices
        )

    with pytest.raises(TypeError):
        mulliken_populations(coeff_ab_mo, occupations, olp_ab_ab, float(num_atoms), ab_atom_indices)

    with pytest.raises(TypeError):
        mulliken_populations(
            coeff_ab_mo, occupations, olp_ab_ab, num_atoms, ab_atom_indices.tolist()
        )
    with pytest.raises(TypeError):
        mulliken_populations(
            coeff_ab_mo, occupations, olp_ab_ab, num_atoms, ab_atom_indices[:, None]
        )
    with pytest.raises(TypeError):
        mulliken_populations(
            coeff_ab_mo, occupations, olp_ab_ab, num_atoms, ab_atom_indices.astype(float)
        )

    with pytest.raises(ValueError):
        mulliken_populations(
            coeff_ab_mo, occupations, olp_ab_ab.reshape(10, 40), num_atoms, ab_atom_indices
        )
    with pytest.raises(ValueError):
        mulliken_populations(
            coeff_ab_mo.reshape(15, 20), occupations, olp_ab_ab, num_atoms, ab_atom_indices
        )
    with pytest.raises(ValueError):
        mulliken_populations(coeff_ab_mo, np.random.rand(20), olp_ab_ab, num_atoms, ab_atom_indices)

    with pytest.raises(TypeError):
        mulliken_populations(
            coeff_ab_mo, occupations, np.random.rand(20, 20), num_atoms, ab_atom_indices.tolist()
        )
    with pytest.raises(TypeError):
        rand_olp_ab_ab = np.random.rand(20, 20)
        rand_olp_ab_ab += rand_olp_ab_ab.T
        mulliken_populations(
            coeff_ab_mo, occupations, rand_olp_ab_ab, num_atoms, ab_atom_indices.tolist()
        )
    with pytest.raises(TypeError):
        rand_olp_ab_ab = np.random.rand(20, 20)
        rand_olp_ab_ab += rand_olp_ab_ab.T
        rand_norm = np.diag(rand_olp_ab_ab) ** (-0.5)
        rand_olp_ab_ab *= rand_norm[:, None]
        rand_olp_ab_ab *= rand_norm[None, :]
        mulliken_populations(
            coeff_ab_mo, occupations, rand_olp_ab_ab, num_atoms, ab_atom_indices.tolist()
        )

    with pytest.raises(ValueError):
        rand_occupations = np.random.rand(15)
        rand_occupations[6] = -1
        mulliken_populations(coeff_ab_mo, rand_occupations, olp_ab_ab, num_atoms, ab_atom_indices)

    # not sure how to check that the warning is raised but the following code prints the warning
    rand_occupations = np.random.rand(15)
    rand_occupations[6] = 2
    mulliken_populations(coeff_ab_mo, rand_occupations, olp_ab_ab, num_atoms, ab_atom_indices)

    with pytest.raises(ValueError):
        rand_weights = np.random.rand(3, 20, 20)
        mulliken_populations(
            coeff_ab_mo,
            rand_occupations,
            olp_ab_ab,
            num_atoms,
            ab_atom_indices,
            atom_weights=rand_weights,
        )
    with pytest.raises(ValueError):
        rand_weights = np.random.rand(4, 20, 19)
        mulliken_populations(
            coeff_ab_mo,
            rand_occupations,
            olp_ab_ab,
            num_atoms,
            ab_atom_indices,
            atom_weights=rand_weights,
        )
    with pytest.raises(ValueError):
        rand_weights = np.random.rand(4, 20, 20)
        mulliken_populations(
            coeff_ab_mo,
            rand_occupations,
            olp_ab_ab,
            num_atoms,
            ab_atom_indices,
            atom_weights=rand_weights,
        )
    with pytest.raises(ValueError):
        rand_weights = np.random.rand(4, 20, 20)
        rand_weights += np.swapaxes(rand_weights, 1, 2)
        mulliken_populations(
            coeff_ab_mo,
            rand_occupations,
            olp_ab_ab,
            num_atoms,
            ab_atom_indices,
            atom_weights=rand_weights,
        )


def test_mulliken_populations():
    """Test orbtools.mulliken.mulliken_populations."""
    # Model system
    coeff_ab_mo = np.identity(10)
    occupations = np.array([2] * 4 + [0] * 6)
    olp_ab_ab = np.identity(10)
    num_atoms = 2
    ab_atom_indices = np.array([0, 0, 1, 1, 0, 0, 1, 1, 0, 1])
    assert np.allclose(
        np.array([4, 4]),
        mulliken_populations(coeff_ab_mo, occupations, olp_ab_ab, num_atoms, ab_atom_indices),
    )
    # H2O RHF/STO-3G
    coeff_ab_mo = np.array(
        [
            [
                9.94099882e-01,
                -2.32889095e-01,
                1.65502866e-08,
                1.00235366e-01,
                6.55422174e-16,
                -1.35631600e-01,
                5.67656304e-08,
            ],
            [
                2.67799213e-02,
                8.31788042e-01,
                -9.03020278e-08,
                -5.23423149e-01,
                -3.01062443e-15,
                9.08581133e-01,
                -4.29452063e-07,
            ],
            [
                3.46630004e-03,
                1.03349385e-01,
                -3.46565859e-01,
                6.48259144e-01,
                3.74502802e-15,
                5.83295647e-01,
                5.82525068e-01,
            ],
            [
                2.72896277e-16,
                -2.21764760e-16,
                4.92555323e-16,
                -6.18875097e-15,
                1.00000000e00,
                1.62509738e-16,
                -6.76860125e-17,
            ],
            [
                2.45105601e-03,
                7.30794097e-02,
                4.90116062e-01,
                4.58390414e-01,
                2.54192740e-15,
                4.12453695e-01,
                -8.23811720e-01,
            ],
            [
                -6.08393842e-03,
                1.60223990e-01,
                4.41542336e-01,
                2.69085788e-01,
                1.34551715e-15,
                -8.07337352e-01,
                8.42614916e-01,
            ],
            [
                -6.08393693e-03,
                1.60223948e-01,
                -4.41542341e-01,
                2.69085849e-01,
                1.99471214e-15,
                -8.07337875e-01,
                -8.42614243e-01,
            ],
        ]
    )
    occupations = np.array([2, 2, 2, 2, 2, 0, 0])
    olp_ab_ab = np.array(
        [
            [1.0, 0.23670392, 0.0, 0.0, 0.0, 0.05490733, 0.05490732],
            [0.23670392, 1.0, 0.0, 0.0, 0.0, 0.47954331, 0.47954323],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.37329955],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.39594342, -0.13197959],
            [0.05490733, 0.47954331, 0.0, 0.0, 0.39594342, 1.0, 0.23846113],
            [0.05490732, 0.47954323, 0.37329955, 0.0, -0.13197959, 0.23846113, 1.0],
        ]
    )
    num_atoms = 3
    ab_atom_indices = np.array([0, 0, 0, 0, 0, 1, 2])
    assert np.allclose(
        np.array([-0.38189777, 0.1909489, 0.19094886]),
        np.array([8, 1, 1])
        - mulliken_populations(coeff_ab_mo, occupations, olp_ab_ab, num_atoms, ab_atom_indices),
    )


def test_mulliken_populations_newbasis():
    """Test orbtools.mulliken.mulliken_populations_newabasis."""
    current_dir = os.path.dirname(__file__)
    coeff_ab_mo = np.load(os.path.join(current_dir, "naclo4_coeff_ab_mo.npy"))
    olp_ab_ab = np.load(os.path.join(current_dir, "naclo4_olp_ab_ab.npy"))
    occupations = np.load(os.path.join(current_dir, "naclo4_occupations.npy"))
    ab_atom_indices = np.load(os.path.join(current_dir, "naclo4_ab_atom_indices.npy"))

    assert np.allclose(
        mulliken_populations_newbasis(
            coeff_ab_mo, occupations, olp_ab_ab, 6, np.identity(124), ab_atom_indices
        ),
        mulliken_populations(coeff_ab_mo, occupations, olp_ab_ab, 6, ab_atom_indices),
    )

    coeff_ab_rand = np.linalg.svd(np.random.rand(124, 124))[0].T
    rand_atom_indices = (np.random.rand(124) * 6 // 1).astype(int)
    # normalize
    olp_rand_rand = coeff_ab_rand.T.dot(olp_ab_ab).dot(coeff_ab_rand)
    coeff_ab_rand *= np.diag(olp_rand_rand) ** (-0.5)
    olp_rand_rand = coeff_ab_rand.T.dot(olp_ab_ab).dot(coeff_ab_rand)
    coeff_rand_mo = project(olp_rand_rand, coeff_ab_rand.T.dot(olp_ab_ab).dot(coeff_ab_mo))
    assert np.allclose(coeff_ab_rand.dot(coeff_rand_mo), coeff_ab_mo)
    assert np.allclose(
        mulliken_populations_newbasis(
            coeff_ab_mo, occupations, olp_ab_ab, 6, coeff_ab_rand, rand_atom_indices
        ),
        mulliken_populations(coeff_rand_mo, occupations, olp_rand_rand, 6, rand_atom_indices),
    )
