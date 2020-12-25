import math


def vdc(k, base=2):
    """[summary]

    Arguments:
        k (int): number

    Keyword Arguments:
        base (int): [description] (default: {2})

    Returns:
        int: [description]
    """
    vdc, denom = 0.0, 1.0
    while k:
        denom *= base
        k, remainder = divmod(k, base)
        vdc += remainder / denom
    return vdc


class vdcorput:
    def __init__(self, base=2):
        self.base = base
        self.count = 0

    def __next__(self):
        """Get the next item

        Returns:
            float:  the next item
        """
        self.count += 1
        return vdc(self.count, self.base)


class halton:
    """Generate base-b Halton sequence

    Arguments:
        n (int): [description]
        b ([int]): sequence base, integer exceeding 1

    Returns:
        ([float]): base-b low discrepancy sequence
    """
    def __init__(self, b):
        self.vdc0 = vdcorput(b[0])
        self.vdc1 = vdcorput(b[1])

    def __next__(self):
        """Get the next item

        Returns:
            list(float):  the next item
        """
        return [next(self.vdc0), next(self.vdc1)]


class circle:
    """Generate Circle Halton sequence 0,..,k

    Arguments:
        k (int): maximum sequence index, non-negative integer

    Keyword Arguments:
        base (int): [description] (default: {2})

    Returns:
        ([float]): base-b low discrepancy sequence
    """
    def __init__(self, base=2):
        self.vc = vdcorput(base)
        self.twopi = 2 * math.pi

    def __next__(self):
        """Get the next item

        Raises:
            StopIteration:  description

        Returns:
            list:  the next item
        """
        vd = next(self.vc)
        theta = self.twopi * vd  # map to [0, 2*math.pi]
        return [math.cos(theta), math.sin(theta)]


class sphere:
    """Generate Sphere Halton sequence 0,..,k

    Arguments:
        k (int): maximum sequence index, non-negative integer

    Keyword Arguments:
        b ([int]): sequence base, integer exceeding 1

    Returns:
        ([float]): base-b low discrepancy sequence
    """
    def __init__(self, b):
        assert len(b) >= 2
        self.cirgen = circle(b[1])
        self.vdc = vdcorput(b[0])

    def __next__(self):
        """Get the next item

        Returns:
            list:  the next item
        """
        vd = next(self.vdc)
        cosphi = 2 * vd - 1  # map to [-1, 1]
        sinphi = math.sqrt(1 - cosphi * cosphi)
        c = next(self.cirgen)
        return [cosphi, sinphi * c[0], sinphi * c[1]]


class sphere3_hopf:
    """
     sphere3_hopf   Halton sequence
     INPUTS   : k - maximum sequence index, non-negative integer
                b - sequence base, integer exceeding 1
    """
    def __init__(self, b):
        assert len(b) >= 3
        self.b = b
        self.vdc0 = vdcorput(b[0])
        self.vdc1 = vdcorput(b[1])
        self.vdc2 = vdcorput(b[2])
        self.twopi = 2 * math.pi

    def __next__(self):
        """Get the next item

        Returns:
            list:  the next item
        """
        vd2 = next(self.vdc2)
        phi = self.twopi * next(self.vdc0)  # map to [0, 2*math.pi]
        psy = self.twopi * next(self.vdc1)  # map to [0, 2*math.pi]
        z = 2 * vd2 - 1  # map to [-1., 1.]
        eta = math.acos(z) / 2
        cos_eta = math.cos(eta)
        sin_eta = math.sin(eta)
        return [
            cos_eta * math.cos(psy), cos_eta * math.sin(psy),
            sin_eta * math.cos(phi + psy),
            sin_eta * math.sin(phi + psy)
        ]

