import numpy as np
from pytest import approx
from scipy.spatial import ConvexHull

from pylds.discrep_2 import discrep_2
from pylds.lds_n import CylinN, SphereN


def rupylds(spgen):
    npoints = 600
    Triples = np.array([spgen.pop() for _ in range(npoints)])
    hull = ConvexHull(Triples)
    triangles = hull.simplices
    return discrep_2(triangles, Triples)


def test_sphere_n():
    spgen = SphereN([2, 3, 5, 7])
    measure = rupylds(spgen)
    assert measure == approx(0.9125914)
    # assert measure < 0.913
    # assert measure > 0.912


def test_cylin_n():
    cygen = CylinN([2, 3, 5, 7])
    measure = rupylds(cygen)
    assert measure == approx(1.0505837105828988)
    # assert measure < 1.086
    # assert measure > 1.085
