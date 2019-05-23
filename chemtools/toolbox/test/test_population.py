"""Test chemtools.toolbox.motbased.MOTBasedTool.compute_populations"""
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path

from chemtools.toolbox.motbased import MOTBasedTool
import numpy as np
from numpy.testing import assert_raises


def test_compute_populations():
    """Test MOTBasedTool.compute_populations.

    Results are compared with those generated from Gaussian. The system is H2O UB3LYP/aug-cc-pVDZ,
    singlet, and at zero net charge. The following the its coordinates:

    O 0.0159484498, 0.0170042791, 0.0238579956
    H -0.772778442, 0.561446550, 1.57501231
    H 1.29850109, 1.26951236, -0.309113326

    """
    with path("chemtools.data.examples", "h2o.fchk") as fname:
        mot = MOTBasedTool.from_file(str(fname))
    assert np.allclose(
        np.array([-0.166945, 0.083473, 0.083473]),
        mot.compute_populations(),
        atol=1e-6,
    )
    assert np.allclose(
        np.array([0.203965, -0.101983, -0.101983]),
        mot.compute_populations("lowdin"),
        atol=1e-6,
    )
    assert_raises(ValueError, mot.compute_populations, "bad type")
