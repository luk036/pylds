from pytest import approx

import numpy as np
from scipy.spatial import ConvexHull

from pylds.discrep_2 import discrep_2
from pylds.low_discr_seq import sphere


def test_sphere():
    npoints = 600
    sgen = sphere([2, 3])
    Triples = np.array([sgen() for _ in range(npoints)])
    hull = ConvexHull(Triples)
    triangles = hull.simplices
    measure = discrep_2(triangles, Triples)
    assert measure == approx(0.2883404521)
    # assert measure < 0.2884
    # assert measure > 0.2883
