import numpy as np
from pytest import approx
from scipy.spatial import ConvexHull

from pylds.discrep_2 import discrep_2
from pylds.low_discr_seq import sphere3_hopf
from pylds.low_discr_seq_n import sphere3


def rupylds(spgen):
    npoints = 600
    Triples = np.array([spgen() for _ in range(npoints)])
    hull = ConvexHull(Triples)
    triangles = hull.simplices
    return discrep_2(triangles, Triples)


def test_sphere():
    spgen = sphere3([2, 3, 5, 7])
    measure = rupylds(spgen)
    assert measure == approx(0.6501446)
    # assert measure < 0.6502
    # assert measure > 0.6501


def test_sphere_hopf():
    spgen = sphere3_hopf([2, 3, 5, 7])
    measure = rupylds(spgen)
    assert measure == approx(0.740866)
    # assert measure < 0.7409
    # assert measure > 0.7408
