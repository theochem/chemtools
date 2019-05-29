"""Test chemtools.utils.mesh."""
from chemtools.utils.mesh import plane_mesh
import numpy as np
from numpy.testing import assert_raises


def test_plane_mesh():
    """Test chemtools.utils.mesh.plane_mesh."""
    assert_raises(TypeError, plane_mesh, np.random.rand(3, 2), 1.0, 1.0)
    assert_raises(TypeError, plane_mesh, np.random.rand(3, 3).tolist(), 1.0, 1.0)
    assert_raises(TypeError, plane_mesh, np.random.rand(3, 3), 1, 1.0)
    assert_raises(TypeError, plane_mesh, np.random.rand(3, 3), 0.0, 1.0)
    assert_raises(TypeError, plane_mesh, np.random.rand(3, 3), -1.0, 1.0)
    assert_raises(TypeError, plane_mesh, np.random.rand(3, 3), 1.0, [])
    assert_raises(TypeError, plane_mesh, np.random.rand(3, 3), 1.0, np.array(1))
    assert_raises(TypeError, plane_mesh, np.random.rand(3, 3), 1.0, -0.1)

    # equilateral triangle on xy plane, centered at origin, side of length 1
    coords = np.array(
        [
            [2.0 / 3 * 0.75 ** 0.5, 0, 0],
            [-1.0 / 3 * 0.75 ** 0.5, 0.5, 0],
            [-1.0 / 3 * 0.75 ** 0.5, -0.5, 0],
        ]
    )
    # make mesh that spans from 1 to - 0.5 in one of the direction, and other side scaled
    # appropriately (3 / sqrt(2))
    ref = []
    height = 1.5
    width = 3 ** 0.5
    for i in range(5):
        temp = []
        for j in range(5):
            temp.append([1 - i * height / 4, -width / 2 + width / 4 * j, 0])
        ref.append(temp)

    assert np.allclose(plane_mesh(coords, 0.5, 1 - 2.0 / 3 * 0.75 ** 0.5), ref)
