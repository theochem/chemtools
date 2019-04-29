import numpy as np
from quasibasis.mulliken import Mulliken

def check_if_exception_raised(func, exception):
    """ Passes if given exception is raised

    Parameter
    ---------
    func : function object
        Function that contains the desired test code
    exception : Exception
        Exception that is raised

    Returns
    -------
    bool
        True if Exception is raised
        False if Exception is not raised
    """
    try:
        func()
    except exception:
        return True
    except:
        return False
    else:
        return False

def test_init():
    coeff_ab_mo = np.random.rand(10, 10)
    occupations = np.array([1]*4+[0]*6)
    ind_occ = np.array([True]*4+[False]*6)
    olp_ab_ab = np.random.rand(10, 10)
    num_center = 2
    basis_map = [0]*2 + [1]*8
    # normalize and symmeterize olp
    olp_ab_ab = olp_ab_ab + olp_ab_ab.T
    norm = (np.diag(olp_ab_ab)**(-0.5)).reshape((1, coeff_ab_mo.shape[0]))
    olp_ab_ab *= norm * norm.T
    # normalize mo
    olp_mo_mo = coeff_ab_mo.T.dot(olp_ab_ab).dot(coeff_ab_mo)
    coeff_ab_mo *= (np.diag(olp_mo_mo)**(-0.5)).reshape((1, coeff_ab_mo.shape[1]))
    mul = Mulliken(coeff_ab_mo, occupations, olp_ab_ab, num_center, basis_map)
    assert np.allclose(mul.coeff_ab_mo, coeff_ab_mo[:, ind_occ])
    assert np.allclose(mul.occupations, occupations[ind_occ])
    assert np.allclose(mul.olp_ab_ab, olp_ab_ab)
    assert mul.num_center == 2
    assert mul.basis_map == tuple(basis_map)
    weight_zero = np.zeros(olp_ab_ab.shape)
    weight_zero[:2] += 0.5
    weight_zero[:, :2] += 0.5
    assert np.allclose(mul.weights[0], weight_zero)
    weight_one = np.zeros(olp_ab_ab.shape)
    weight_one[2:] += 0.5
    weight_one[:, 2:] += 0.5
    assert np.allclose(mul.weights[1], weight_one)
   # with weights
    weights = [np.ones(olp_ab_ab.shape)*0.5]*2
    mul = Mulliken(coeff_ab_mo, occupations, olp_ab_ab, num_center, basis_map, weights=weights)
    assert np.allclose(mul.weights[1], np.zeros(olp_ab_ab.shape)+0.5)
    assert np.allclose(mul.weights[1], np.zeros(olp_ab_ab.shape)+0.5)
    # property
    assert mul.num_ab == coeff_ab_mo.shape[0] == olp_ab_ab.shape[0] == olp_ab_ab.shape[1]

def test_weights_default():
    coeff_ab_mo = np.random.rand(10, 10)
    occupations = np.array([1]*4+[0]*6)
    olp_ab_ab = np.random.rand(10, 10)
    num_center = 2
    basis_map = [0]*2 + [1]*8
    # normalize and symmeterize olp
    olp_ab_ab = olp_ab_ab + olp_ab_ab.T
    norm = (np.diag(olp_ab_ab)**(-0.5)).reshape((1, coeff_ab_mo.shape[0]))
    olp_ab_ab *= norm * norm.T
    # normalize mo
    olp_mo_mo = coeff_ab_mo.T.dot(olp_ab_ab).dot(coeff_ab_mo)
    coeff_ab_mo *= (np.diag(olp_mo_mo)**(-0.5)).reshape((1, coeff_ab_mo.shape[1]))
    mul = Mulliken(coeff_ab_mo, occupations, olp_ab_ab, num_center, basis_map)
    # test
    weights = mul.weights_default(num_center, basis_map)
    weight_zero = np.zeros(olp_ab_ab.shape)
    weight_zero[:2] += 0.5
    weight_zero[:, :2] += 0.5
    assert np.allclose(weights[0], weight_zero)
    weight_one = np.zeros(olp_ab_ab.shape)
    weight_one[2:] += 0.5
    weight_one[:, 2:] += 0.5
    assert np.allclose(weights[1], weight_one)

def test_get_population():
    """ Tests Mulliken.get_population

    """
    # Model system
    coeff_ab_mo = np.identity(10)
    occupations = np.array([2]*4 + [0]*6)
    olp_ab_ab = np.identity(10)
    num_center = 2
    basis_map = [0, 0, 1, 1, 0, 0, 1, 1, 0, 1]
    mul = Mulliken(coeff_ab_mo, occupations, olp_ab_ab, num_center, basis_map)
    assert np.allclose(np.array([4, 4]), mul.get_population())
    # H2O RHF/STO-3G
    coeff_ab_mo = np.array([[ 9.94099882e-01,-2.32889095e-01, 1.65502866e-08,
                              1.00235366e-01, 6.55422174e-16,-1.35631600e-01,
                              5.67656304e-08],
                            [ 2.67799213e-02, 8.31788042e-01,-9.03020278e-08,
                             -5.23423149e-01,-3.01062443e-15, 9.08581133e-01,
                             -4.29452063e-07],
                            [ 3.46630004e-03, 1.03349385e-01,-3.46565859e-01,
                              6.48259144e-01, 3.74502802e-15, 5.83295647e-01,
                              5.82525068e-01],
                            [ 2.72896277e-16,-2.21764760e-16, 4.92555323e-16,
                             -6.18875097e-15, 1.00000000e+00, 1.62509738e-16,
                             -6.76860125e-17],
                            [ 2.45105601e-03, 7.30794097e-02, 4.90116062e-01,
                              4.58390414e-01, 2.54192740e-15, 4.12453695e-01,
                             -8.23811720e-01],
                            [-6.08393842e-03, 1.60223990e-01, 4.41542336e-01,
                              2.69085788e-01, 1.34551715e-15,-8.07337352e-01,
                              8.42614916e-01],
                            [-6.08393693e-03, 1.60223948e-01,-4.41542341e-01,
                              2.69085849e-01, 1.99471214e-15,-8.07337875e-01,
                             -8.42614243e-01]])
    occupations = np.array([2, 2, 2, 2, 2, 0, 0])
    olp_ab_ab = np.array([[ 1., 0.23670392, 0., 0., 0., 0.05490733, 0.05490732],
                        [ 0.23670392, 1., 0., 0., 0., 0.47954331, 0.47954323],
                        [ 0., 0., 1., 0., 0., 0., 0.37329955],
                        [ 0., 0., 0., 1., 0., 0. , 0. ],
                        [ 0., 0., 0., 0., 1., 0.39594342, -0.13197959],
                        [ 0.05490733, 0.47954331, 0., 0., 0.39594342, 1., 0.23846113],
                        [ 0.05490732, 0.47954323, 0.37329955, 0., -0.13197959, 0.23846113, 1.]])
    num_center = 3
    basis_map = [0, 0, 0, 0, 0, 1, 2]
    mul = Mulliken(coeff_ab_mo, occupations, olp_ab_ab, num_center, basis_map)
    pop = mul.get_population()
    partial_pop = np.array([8, 1, 1]) - pop
    assert np.allclose(np.array([-0.38189777,  0.1909489 ,  0.19094886]),
                       partial_pop)


test_init()
test_weights_default()
test_get_population()
