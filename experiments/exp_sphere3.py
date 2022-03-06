from __future__ import print_function

from pprint import pprint

import matplotlib.pylab as lab
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull

from pylds.discrep_2 import discrep_2
from pylds.low_discr_seq import sphere3_hopf
from pylds.low_discr_seq_n import sphere3


def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec.transpose()


def dispersion(Triples):
    hull = ConvexHull(Triples)
    triangles = hull.simplices
    measure = discrep_2(triangles, Triples)
    return measure


npoints = 4000
ndim = 4
Triples_r = sample_spherical(npoints, ndim)
# Triples_h = np.array([p for p in sphere3_hopf(npoints, [2, 3, 5])])
sh = sphere3_hopf([2, 3, 5])
Triples_h = np.array([sh() for _ in range(npoints)])
s3 = sphere3([2, 3, 5])
# Triples_s = np.array([p for p in sphere3(npoints, [2, 3, 5])])
Triples_s = np.array([s3() for _ in range(npoints)])

x = list(range(100, npoints, 100))
res_r = []
res_h = []
res_s = []

for i in x:
    res_r += [dispersion(Triples_r[:i, :])]
    res_h += [dispersion(Triples_h[:i, :])]
    res_s += [dispersion(Triples_s[:i, :])]

plt.plot(x, res_r, "r", label="Random")
plt.plot(x, res_h, "b", label="Hopf")
plt.plot(x, res_s, "g", label="Our")
plt.legend(loc="best")
plt.xlabel("#points")
plt.ylabel("discrepancy")
plt.show()
