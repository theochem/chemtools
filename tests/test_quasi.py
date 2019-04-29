import numpy as np
from quasibasis.quasi import QuasiTransformation, project
from quasibasis.mulliken import Mulliken
from quasibasis.wrapper_horton import HortonData


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


def normalize(olp, coeff):
    norm = np.diag(coeff.T.dot(olp).dot(coeff))
    return coeff * norm ** (-0.5)


def test_project():
    """ Tests the quasi.project

    """
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


def test_init():
    """ Tests the quasi.QuasiTransformation.__init__ and the various properties

    """
    coeff_prim1_ab = normalize(np.identity(10), np.random.rand(10, 10))
    olp_ab_ab = coeff_prim1_ab.T.dot(np.identity(10)).dot(coeff_prim1_ab)
    coeff_prim2_aao = normalize(np.identity(5), np.random.rand(5, 5))
    olp_aao_aao = coeff_prim2_aao.T.dot(np.identity(5)).dot(coeff_prim2_aao)
    olp_aao_ab = np.random.rand(5, 10)
    coeff_ab_mo = normalize(olp_ab_ab, np.random.rand(10, 10))
    indices_span = np.array([True] * 2 + [False] * 8)
    # initialize
    quasi = QuasiTransformation(coeff_ab_mo, olp_ab_ab, olp_aao_ab, olp_aao_aao, indices_span)
    assert np.allclose(quasi._coeff_ab_mo, coeff_ab_mo)
    assert np.allclose(quasi._olp_ab_ab, olp_ab_ab)
    assert np.allclose(quasi._olp_aao_ab, olp_aao_ab)
    assert np.allclose(quasi._indices_span, indices_span)
    # indices_span setter
    quasi.indices_span = np.array([2, 3])
    assert np.allclose(quasi._indices_span, np.array([False] * 2 + [True] * 2 + [False] * 6))
    quasi.indices_span = np.array([False] * 2 + [True] * 2 + [False] * 6)
    assert np.allclose(quasi._indices_span, np.array([False] * 2 + [True] * 2 + [False] * 6))
    # properties
    assert np.allclose(quasi.coeff_ab_mo, coeff_ab_mo)
    assert np.allclose(quasi.olp_ab_ab, olp_ab_ab)
    assert np.allclose(quasi.olp_aao_ab, olp_aao_ab)
    assert np.allclose(quasi.olp_aao_aao, olp_aao_aao)
    assert np.allclose(quasi.indices_span, np.array([False] * 2 + [True] * 2 + [False] * 6))
    assert np.allclose(quasi.olp_ab_aao, olp_aao_ab.T)
    olp_mo_mo = coeff_ab_mo.T.dot(olp_ab_ab).dot(coeff_ab_mo)
    assert np.allclose(quasi.olp_mo_mo, olp_mo_mo)
    assert np.allclose(quasi.olp_aao_mo, olp_aao_ab.dot(coeff_ab_mo))
    assert np.allclose(quasi.olp_mo_aao, coeff_ab_mo.T.dot(olp_aao_ab.T))
    omo_ind = quasi.indices_span
    assert np.allclose(quasi.coeff_ab_omo, coeff_ab_mo[:, omo_ind])
    assert np.allclose(quasi.olp_omo_omo, olp_mo_mo[omo_ind, :][:, omo_ind])
    assert np.allclose(quasi.olp_omo_ab, (coeff_ab_mo[:, omo_ind]).T.dot(olp_ab_ab))
    assert np.allclose(quasi.olp_ab_omo, olp_ab_ab.dot(coeff_ab_mo[:, omo_ind]))
    assert np.allclose(quasi.olp_omo_aao, (coeff_ab_mo[:, omo_ind]).T.dot(olp_aao_ab.T))
    assert np.allclose(quasi.olp_aao_omo, olp_aao_ab.dot(coeff_ab_mo[:, omo_ind]))
    assert quasi.num_ab == olp_ab_ab.shape[0] == olp_aao_ab.shape[1] == coeff_ab_mo.shape[0]
    assert quasi.num_aao == olp_aao_ab.shape[0] == olp_aao_aao.shape[0]
    assert quasi.num_mo == coeff_ab_mo.shape[1]


def test_quambo():
    """ Tests QUasiTransformation.quambo by comparing to literature value
    (population analysis with QUAMBO basis)
    """
    hd = HortonData("naclo4_631gdp_hf_g09.fchk", "minao.gbs")
    coeff_ab_mo_sep = hd.coeff_ab_mo_sep
    olp_ab_ab_sep = hd.olp_ab_ab_sep
    olp_cao_ab_sep = hd.olp_cao_ab_sep
    olp_cao_cao_sep = hd.olp_cao_cao_sep
    occupations_sep = hd.occupations_sep
    basis_map_sep = hd.cao_basis_map_sep
    quasi = QuasiTransformation(
        coeff_ab_mo_sep[0],
        olp_ab_ab_sep[0],
        olp_cao_ab_sep[0],
        olp_cao_cao_sep[0],
        occupations_sep[0].astype(bool),
    )
    coeff_ab_quambo = quasi.quambo()
    olp_quambo_quambo = coeff_ab_quambo.T.dot(quasi.olp_ab_ab).dot(coeff_ab_quambo)
    olp_quambo_omo = coeff_ab_quambo.T.dot(quasi.olp_ab_omo)
    coeff_quambo_omo = project(olp_quambo_quambo, olp_quambo_omo)
    occupations = np.array([i for i in occupations_sep[0] if i > 0])
    pop = Mulliken(
        coeff_quambo_omo, occupations, olp_quambo_quambo, 6, basis_map_sep[0]
    ).get_population()
    partial_pop = np.array([11, 17, 8, 8, 8, 8]) - pop
    assert np.allclose(
        partial_pop, np.array([0.968, 2.454, -0.807, -0.904, -0.904, -0.807]), atol=1e-3
    )


def test_quao():
    """ Tests QUasiTransformation.quao by comparing to literature value
    (population analysis with QUAO basis)
    """
    hd = HortonData("naclo4_631gdp_hf_g09.fchk", "minao.gbs")
    coeff_ab_mo_sep = hd.coeff_ab_mo_sep
    olp_ab_ab_sep = hd.olp_ab_ab_sep
    olp_cao_ab_sep = hd.olp_cao_ab_sep
    olp_cao_cao_sep = hd.olp_cao_cao_sep
    occupations_sep = hd.occupations_sep
    basis_map_sep = hd.cao_basis_map_sep
    quasi = QuasiTransformation(
        coeff_ab_mo_sep[0],
        olp_ab_ab_sep[0],
        olp_cao_ab_sep[0],
        olp_cao_cao_sep[0],
        occupations_sep[0].astype(bool),
    )
    coeff_ab_quao = quasi.quao()
    olp_quao_quao = coeff_ab_quao.T.dot(quasi.olp_ab_ab).dot(coeff_ab_quao)
    olp_quao_omo = coeff_ab_quao.T.dot(quasi.olp_ab_omo)
    coeff_quao_omo = project(olp_quao_quao, olp_quao_omo)
    occupations = np.array([i for i in occupations_sep[0] if i > 0])
    pop = Mulliken(coeff_quao_omo, occupations, olp_quao_quao, 6, basis_map_sep[0]).get_population()
    partial_pop = np.array([11, 17, 8, 8, 8, 8]) - pop
    assert np.allclose(
        partial_pop, np.array([0.967, 2.498, -0.819, -0.914, -0.914, -0.819]), atol=1e-3
    )


def test_iao():
    """ Tests QUasiTransformation.iao by comparing to literature value
    (population analysis with IAO basis)
    """
    hd = HortonData("ch4_svp_minao_iao.fchk", "minao.gbs")
    coeff_ab_mo_sep = hd.coeff_ab_mo_sep
    olp_ab_ab_sep = hd.olp_ab_ab_sep
    olp_cao_ab_sep = hd.olp_cao_ab_sep
    olp_cao_cao_sep = hd.olp_cao_cao_sep
    occupations_sep = hd.occupations_sep
    basis_map_sep = hd.cao_basis_map_sep
    quasi = QuasiTransformation(
        coeff_ab_mo_sep[0],
        olp_ab_ab_sep[0],
        olp_cao_ab_sep[0],
        olp_cao_cao_sep[0],
        occupations_sep[0].astype(bool),
    )
    coeff_ab_iao = quasi.iao()
    olp_iao_iao = coeff_ab_iao.T.dot(quasi.olp_ab_ab).dot(coeff_ab_iao)
    olp_iao_omo = coeff_ab_iao.T.dot(quasi.olp_ab_omo)
    coeff_iao_omo = project(olp_iao_iao, olp_iao_omo)
    occupations = np.array([i for i in occupations_sep[0] if i > 0])
    pop = Mulliken(coeff_iao_omo, occupations, olp_iao_iao, 5, basis_map_sep[0]).get_population()
    partial_pop = np.array([6, 1, 1, 1, 1]) - pop
    assert np.allclose(partial_pop, np.array([-0.49, 0.12, 0.12, 0.12, 0.12]), atol=1e-2)


test_project()
test_init()
test_quambo()
test_quao()
test_iao()
