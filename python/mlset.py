"""
Convenience layer module for integration using multiple level sets.
"""
from xfem import *
from itertools import chain, permutations, product


class DomainTypeArray():
    """
    Container class for multiple level sets.

    Parameters
    ----------
    dtlist : list(tuples) or singe tuple
        The tuple or a list of tuples of DOMAIN_TYPE describing the 
        region where we integrate. ANY is a not valid DOMAIN_TYPE but 
        is expanded on initialisation onto POS and NEG.


    Attributes
    ----------
    dtlist : list(tuples)
        Expanded list containing valid tuples describing the region of
        interest
    codim : int
        Co-dimension of the  region of interest.


    Methods
    -------
    Boundary(element_marking=False):
        - A DomainTypeArray describing the boundary of the current 
        DomainTypeArray region.


    """

    def __init__(self, dtlist):
        if type(dtlist) == tuple:
            dtlist = [dtlist]

        if type(dtlist) != list or type(dtlist[0]) != tuple:
            raise TypeError("Invalid input: dtlist must be list of tuples")

        if (ANY in chain.from_iterable(dtlist)):
            dtlist_expand = []
            for dtt in dtlist:
                if ANY in dtt:
                    locations_expand = [i for i in range(len(dtt)) if dtt[i] == ANY]
                    n = len(locations_expand)
                    dt_lists = [list(dtt) for i in range(2**n)]
                    dt_per = list(set(permutations([NEG] * n + [POS] * n, n)))
                    for i, dtl in enumerate(dt_lists):
                        for j, k in enumerate(locations_expand):
                            dt_lists[i][k] = dt_per[i][j]
                        dtlist_expand.append(tuple(dt_lists[i]))
                else:
                    dtlist_expand.append(dtt)
            self.dtlist = list(set(dtlist_expand))
        else:
            self.dtlist = list(set(dtlist))

        self.codim = 0
        for dt in self.dtlist[0]:
            if dt == IF:
                self.codim += 1

        # Safety checks
        if self.codim > 3:
            print("Warning: Codim > 3 !!!")

        dtt_len = len(self.dtlist[0])
        for dtt in self.dtlist:
            for dt in dtt:
                if dt not in [NEG, POS, IF, ANY]:
                    raise Exception("Initialised with something other than NEG, POS, IF, ANY")
            n = dtt.count(IF)
            if n != self.codim:
                raise Exception("DomainTypeArray initialised with tuples of "
                                "different co-dimension!")
            if dtt_len != len(dtt):
                raise Exception("DomainTypeArray initialised with tuples of "
                                "different length!")

    def __len__(self):
        return self.dtlist.__len__()

    def __iter__(self):
        return self.dtlist.__iter__()

    def __contains__(self, dtt):
        return dtt in self.dtlist

    def __str__(self):
        str_out = "DomainTypeArray:\n"
        str_out += "  Co-dim : {}\n".format(self.codim)
        str_out += "  Domains: ["
        for dtt in self.dtlist:
            str_out += dtt.__str__() + ", "
        return str_out[:-2] + "]"

    def __or__(self, dta_b):
        if self.codim != dta_b.codim:
            raise Exception("Co-dims don't match: Union not possible!")
        return DomainTypeArray(self.dtlist + dta_b.dtlist)

    def __and__(self, dta_b):
        if self.codim == 0:
            return DomainTypeArray([dtt for dtt in self.dtlist if dtt in dta_b.dtlist])
        else:
            dtl_out = []
            for dtt1, dtt2 in product(self.dtlist, dta_b.dtlist):
                if dtt1 == dtt2:
                    dtl_out.append(dtt1)
                    continue
                dtt_out = []
                for i, dt in enumerate(dtt1):
                    if dt != dtt2[i]:
                        if dt in [POS, NEG] and dtt2[i] in [POS, NEG]:
                            break
                        else:
                            dtt_out.append(IF)
                    else:
                        dtt_out.append(dt)
                else:
                    dtl_out.append(tuple(dtt_out))
            return DomainTypeArray(dtl_out)

    def __invert__(self):

        n = len(self.dtlist[0])
        dt_per = permutations([NEG] * (n - self.codim) +
                              [POS] * (n - self.codim) + [IF] * self.codim, n)
        dt_per = list(filter(lambda dtt: dtt.count(IF) == self.codim,
                             set(dt_per)))
        for dt in self.dtlist:
            dt_per.remove(dt)
        return DomainTypeArray(dt_per)

    def __ior__(self, dta_b):
        return self.__or__(dta_b)

    def __iand__(self, dta_b):
        return self.__and__(dta_b)

    def Boundary(self):
        """
        Compute the Boundary of the current region.

        Returns
        -------
        DomainTypeArray
        """
        if self.codim >= 3:
            raise Exception("Boundary does not make sense for codim >=3")

        dtl_out = []
        for dtt in self.dtlist:
            for i, dt in enumerate(dtt):
                dtt_out = list(dtt)
                if dt != IF:
                    dtt_out[i] = IF
                    dtl_out.append(tuple(dtt_out))
        return DomainTypeArray(dtl_out)


def dta_indicator(lsets, dta):
    """
    Indicator function for a DomainTypeArray

    Parameters
    ----------
    lsets : tuple(CoefficientFunctions)
        The level set functions defining the region
    dta : DomainTypeArray
        The dta describing the region of interest. Must have codim=0

    Returns
    -------
    CoefficientFunction
        1 in the region on interest, 0 else
    """

    if dta.codim != 0:
        raise NotImplementedError("Indicator only for codim=1")

    ind = CoefficientFunction(1)
    for dtt in dta.dtlist:
        for i, dt in enumerate(dtt):
            if dt == POS:
                ind *= IfPos(lsets[i], CoefficientFunction(1),
                             CoefficientFunction(0))
            elif dt == NEG:
                ind *= IfPos(-lsets[i], CoefficientFunction(1),
                             CoefficientFunction(0))

    return ind


if __name__ == "__main__":

    dta1 = DomainTypeArray([(POS, NEG)])
    dta2 = DomainTypeArray([(POS, POS)])
    dta3 = dta1 | dta2

    print((POS, NEG) in dta3)
    print(len(dta3))
    input("")

    print(dta3)
    print("Expected: ", [(POS, NEG), (POS, POS)])
    input("")

    dta3 |= DomainTypeArray([(NEG, POS)])

    print(dta3)
    print("Expected: ", [(POS, NEG), (POS, POS), (NEG, POS)])
    input("")

    dta4 = dta3 & dta1

    print(dta4)
    print("Expected: ", [(POS, NEG)])
    input("")

    dta3 &= dta1 | dta2
    print(dta3)
    print("Expected: ", [(POS, NEG), (POS, POS)])
    input("")

    dta5 = DomainTypeArray([(IF, NEG)]) & DomainTypeArray([(POS, NEG)])
    print(dta5)
    print("Expected: ", [(IF, NEG)])
    input("")

    dta5 &= DomainTypeArray([(POS, IF)])
    print(dta5)
    print("Expected: ", [(IF, IF)])
    input("")

    print(~dta3)
    print("Expected: ", [(NEG, POS), (NEG, NEG)])
    input("")

    dta6 = DomainTypeArray([(IF, NEG, POS)])
    dta6 &= DomainTypeArray([(POS, IF, POS)])
    print(~dta6)
    print("Expected: ", [(IF, NEG, IF), (POS, IF, IF), (IF, POS, IF),
                         (NEG, IF, IF), (IF, IF, NEG)])
    input("")

    dta7 = DomainTypeArray([(NEG, NEG, NEG)])
    print("Boundary: ", dta7.Boundary())
    print("Expected: ", [(NEG, IF, NEG), (NEG, NEG, IF), (IF, NEG, NEG)])
    input("")

    print("elements_makings=True\nBoundary: ",
          dta7.Boundary(element_marking=True))
    print("Expected: ", [(NEG, IF, NEG), (NEG, NEG, IF), (IF, NEG, NEG),
                         (IF, IF, NEG), (IF, NEG, IF), (NEG, IF, IF),
                         (IF, IF, IF)])
    input("")

    dta8 = DomainTypeArray([(POS, IF, NEG)])
    print("Boundary: ", dta8.Boundary())
    print("Expected: ", [(POS, IF, IF), (IF, IF, NEG)])
    input("")
