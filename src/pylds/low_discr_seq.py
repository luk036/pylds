import math
from typing import List

twoPI = 2 * math.pi


def vdc(k: int, base: int = 2) -> float:
    """[summary]

    Arguments:
        k (int): number

    Keyword Arguments:
        base (int): [description] (default: {2})

    Returns:
        int: [description]
    """
    vdc = 0.0
    denom = 1.0
    while k != 0:
        denom *= base
        remainder = k % base
        k //= base
        vdc += remainder / denom
    return vdc


class vdcorput:
    def __init__(self, base: int = 2):
        """[summary]

        Args:
            base (int, optional): [description]. Defaults to 2.
        """
        self.base: int = base
        self.count: int = 0

    def __call__(self) -> float:
        """[summary]

        Returns:
            float: [description]
        """
        self.count += 1
        return vdc(self.count, self.base)

    def reseed(self, seed: int):
        self.count = seed


class halton:
    """Generate Halton sequence"""

    def __init__(self, base: List[int]):
        self.vdc0 = vdcorput(base[0])
        self.vdc1 = vdcorput(base[1])

    def __call__(self) -> List[float]:
        """Get the next item

        Returns:
            list(float):  the next item
        """
        return [self.vdc0(), self.vdc1()]

    def reseed(self, seed: int):
        self.vdc0.reseed(seed)
        self.vdc1.reseed(seed)


class circle:
    """Generate Circle Halton sequence 0,..,k

    Arguments:
        k (int): maximum sequence index, non-negative integer

    Keyword Arguments:
        base (int): [description] (default: {2})

    Returns:
        ([float]): base-b low discrepancy sequence
    """

    def __init__(self, base: int = 2):
        self.vdc = vdcorput(base)

    def __call__(self) -> List[float]:
        """Get the next item

        Raises:
            StopIteration:  description

        Returns:
            list:  the next item
        """
        theta = twoPI * self.vdc()  # map to [0, 2*math.pi]
        return [math.cos(theta), math.sin(theta)]

    def reseed(self, seed: int):
        self.vdc.reseed(seed)


class sphere:
    """Generate Sphere Halton sequence 0,..,k

    Arguments:
        k (int): maximum sequence index, non-negative integer

    Keyword Arguments:
        b ([int]): sequence base, integer exceeding 1

    Returns:
        ([float]): base-b low discrepancy sequence
    """

    def __init__(self, base: List[int]):
        assert len(base) >= 2
        self.vdc = vdcorput(base[0])
        self.cirgen = circle(base[1])

    def __call__(self) -> List[float]:
        """Get the next item

        Returns:
            list:  the next item
        """
        cosphi = 2 * self.vdc() - 1  # map to [-1, 1]
        sinphi = math.sqrt(1 - cosphi * cosphi)
        cc = self.cirgen()
        return [cosphi, sinphi * cc[0], sinphi * cc[1]]

    def reseed(self, seed: int):
        self.cirgen.reseed(seed)
        self.vdc.reseed(seed)


class sphere3_hopf:
    """
    sphere3_hopf   Halton sequence
    INPUTS   : k - maximum sequence index, non-negative integer
               b - sequence base, integer exceeding 1
    """

    def __init__(self, base: List[int]):
        assert len(base) >= 3
        self.vdc0 = vdcorput(base[0])
        self.vdc1 = vdcorput(base[1])
        self.vdc2 = vdcorput(base[2])

    def __call__(self) -> List[float]:
        """Get the next item

        Returns:
            list:  the next item
        """
        phi = self.vdc0() * twoPI  # map to [0, 2*math.pi]
        psy = self.vdc1() * twoPI  # map to [0, 2*math.pi]
        # zzz = self.vdc2() * 2 - 1  # map to [-1., 1.]
        # eta = math.acos(zzz) / 2
        # cos_eta = math.cos(eta)
        # sin_eta = math.sin(eta)
        vd = self.vdc2()
        cos_eta = math.sqrt(vd)
        sin_eta = math.sqrt(1 - vd)
        return [
            cos_eta * math.cos(psy),
            cos_eta * math.sin(psy),
            sin_eta * math.cos(phi + psy),
            sin_eta * math.sin(phi + psy),
        ]

    def reseed(self, seed: int):
        self.vdc0.reseed(seed)
        self.vdc1.reseed(seed)
        self.vdc2.reseed(seed)


class halton_n:
    """Generate base-b Halton sequence

    Arguments:
        n (int): [description]
        b ([int]): sequence base, integer exceeding 1

    Returns:
        ([float]): base-b low discrepancy sequence
    """

    def __init__(self, n: int, base: List[int]):
        self.vec_vdc = [vdcorput(base[i]) for i in range(n)]

    def __call__(self) -> List[float]:
        """Get the next item

        Returns:
            list(float):  the next item
        """
        return [vdc() for vdc in self.vec_vdc]

    def reseed(self, seed: int):
        for vdc in self.vec_vdc:
            vdc.reseed(seed)


if __name__ == "__main__":
    halgen = halton_n(4, [2, 5, 7, 3])
    for _ in range(10):
        print(halgen())
