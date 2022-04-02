import numpy as np
from pytest import approx
from scipy.spatial import ConvexHull

from pylds.discrep_2 import discrep_2
from pylds.lds_n import CylinN, SphereN


def run_pylds(spgen):
    npoints = 600
    Triples = np.array([spgen.pop() for _ in range(npoints)])
    hull = ConvexHull(Triples)
    triangles = hull.simplices
    return discrep_2(triangles, Triples)


def test_sphere_n():
    spgen = SphereN(3)
    measure = run_pylds(spgen)
    assert measure == approx(0.9125914)
    # assert measure < 0.913
    # assert measure > 0.912


def test_cylin_n():
    cygen = CylinN(3)
    measure = run_pylds(cygen)
    assert measure == approx(1.0505837105828988)
    # assert measure < 1.086
    # assert measure > 1.085
