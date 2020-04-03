"""
Convenience layer module for integration using multiple level sets.
"""
from ngsolve import Norm, Grad, GridFunction
from xfem import *
from itertools import chain, permutations, product
from collections import Counter


class DomainTypeArray():
    """
    Container class for multiple level sets.

    Parameters
    ----------
    dtlist : list(tuples) or singe tuple
        The tuple or a list of tuples of DOMAIN_TYPE describing the 
        region where we integrate. ANY is a not valid DOMAIN_TYPE but 
        is expanded on initialisation onto POS and NEG.
    lsets : tuple(ngsolve.GridFunction)
        Optional argument: A tuple of discrete level set functions, 
        which are used to compress the instance to only contain regions
        with positive measure with respect to these level sets. If this
        is given, then the DomainTypeArray is compressed on 
        initialisation.

    Attributes
    ----------
    as_list : list(tuples)
        Expanded list containing valid tuples describing the region of
        interest.
    codim : int
        Co-dimension of the  region of interest.
    compressed : bool
        Indicates, whether the instance has been compressed.

    Methods
    -------
    Boundary(element_marking=False):
        - A DomainTypeArray describing the boundary of the current 
        DomainTypeArray region.
    Indicator(lsets):
        Returns an indicator CoefficientFunction of the current 
        (codim=0) DomainTypeArray region.
    IndicatorSmoothed(lsets, eps):
        Returns an indicator CoefficientFunction of an eps-region 
        around the current (codim>0) region.
    GetOuterNormals(lsetsp1):
        Return a dictionary, which contains a dictionary for each region
        which in turn contains the outward pointing unit normal vector
        on each boundary segment. Only for codim=0 DomainTypeArrays.

    """

    def __init__(self, dtlist, lsets=None):
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
            self.as_list = list(set(dtlist_expand))
        else:
            self.as_list = list(set(dtlist))

        self.codim = 0
        for dt in self.as_list[0]:
            if dt == IF:
                self.codim += 1
        self.lsets = lsets
        self.compressed = False
       
        # Safety checks
        if self.codim > 3:
            print("Warning: Codim > 3 !!!")

        dtt_len = len(self.as_list[0])
        for dtt in self.as_list:
            for dt in dtt:
                if dt not in [NEG, POS, IF, ANY]:
                    raise Exception("Initialised with something other than NEG"
                                    ", POS, IF, ANY")
            n = dtt.count(IF)
            if n != self.codim:
                raise Exception("DomainTypeArray initialised with tuples of "
                                "different co-dimension!")
            if dtt_len != len(dtt):
                raise Exception("DomainTypeArray initialised with tuples of "
                                "different length!")
        if self.lsets:
            if len(self.lsets) != dtt_len:
                raise Exception("The number of level sets does not match the "
                    " length of the arrays domain tuples!")
            if type(self.lsets[0]) != GridFunction:
                raise Exception("The level set functions need to be "
                                "ngsolve.GridFunctions!")

        # Compression if possible
        if self.lsets:
            self.Compress(self.lsets)
            self.compressed = True


    def __len__(self):
        return self.as_list.__len__()

    def __iter__(self):
        return self.as_list.__iter__()

    def __contains__(self, dtt):
        return dtt in self.as_list

    def __eq__(self, dta_b):
        if self.codim != dta_b.codim or self.__len__() != len(dta_b):
            return False
        else:
            return Counter(self.as_list) == Counter(dta_b.as_list)

    def __ne__(self, dta_b):
        return not self.__eq__(dta_b)

    def __str__(self):
        str_out = "DomainTypeArray:\n"
        str_out += "  Co-dim : {}\n".format(self.codim)
        str_out += "  Domains: ["
        for dtt in self.as_list:
            str_out += dtt.__str__() + ", "
        return str_out[:-2] + "]"

    def __or__(self, dta_b):
        if self.codim != dta_b.codim:
            raise Exception("Co-dims don't match: Union not possible!")
        return DomainTypeArray(self.as_list + dta_b.as_list, self.lsets)

    def __and__(self, dta_b):
        if self.codim == 0:
            return DomainTypeArray([dtt for dtt in self.as_list if dtt in dta_b.as_list], self.lsets)
        else:
            dtl_out = []
            for dtt1, dtt2 in product(self.as_list, dta_b.as_list):
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
            return DomainTypeArray(dtl_out, self.lsets)

    def __invert__(self):
        n = len(self.as_list[0])
        dt_per = permutations([NEG] * (n - self.codim) +
                              [POS] * (n - self.codim) + [IF] * self.codim, n)
        dt_per = list(filter(lambda dtt: dtt.count(IF) == self.codim,
                             set(dt_per)))
        for dt in self.as_list:
            dt_per.remove(dt)
        return DomainTypeArray(dt_per, self.lsets)

    def __ior__(self, dta_b):
        return self.__or__(dta_b)

    def __iand__(self, dta_b):
        return self.__and__(dta_b)

    def Compress(self, lsets):
        """
        For a given set of level sets, this function removes any regions
        contained in the instance, which do not have a positive measure.
        The level sets remain attached to the instance, such that 
        self.Boundary() called after self.Compress(lsets) automatically
        compresses the output. 
        
        """

        self.lsets = lsets
        self.compressed = True

        if len(lsets) != len(self.as_list[0]):
            raise Exception("The number of level sets does not match the "
                            " length of the arrays domain tuples!")
        if type(self.lsets[0]) != GridFunction:
            raise Exception("The level set functions need to be "
                            "ngsolve.GridFunctions!")

        for dtt in self.as_list:
            weight = weight = Integrate({"levelset": self.lsets, 
                                         "domain_type": dtt},
                                         cf=1, 
                                         mesh=self.lsets[0].space.mesh,
                                         order=0)
            if abs(weight) < 1e-12:
                self.as_list.remove(dtt)


    def Boundary(self):
        """
        Compute the Boundary of the current region. If the instance has
        been compressed, then only level set boundaries with positive 
        measure are returned.

        Returns
        -------
        DomainTypeArray
        """
        if self.codim >= 3:
            raise Exception("Boundary does not make sense for codim >=3")

        dtl_out = []
        for dtt in self.as_list:
            for i, dt in enumerate(dtt):
                if dt != IF:
                    if ( dtt[:i] + tuple([POS]) + dtt[i+1:] in self and 
                         dtt[:i] + tuple([NEG]) + dtt[i+1:] in self):
                        continue
                    else:
                        is_relevant = True
                        dtt_out = dtt[:i] + tuple([IF]) + dtt[i+1:]

                        if self.lsets:
                            weight = Integrate({"levelset": self.lsets,
                                                "domain_type": dtt_out},
                                                cf=1, 
                                                mesh=self.lsets[0].space.mesh,
                                                order=0)
                            if abs(weight) < 1e-12:
                                is_relevant = False
                        if is_relevant:
                            dtl_out.append(dtt_out)
    
        return DomainTypeArray(dtl_out, self.lsets)

    def Indicator(self, lsets):
        """
        Indicator function for a DomainTypeArray of codim=0

        Parameters
        ----------
        lsets : tuple(ngsolve.CoefficientFunctions)
            The level set functions defining the region

        Returns
        -------
        CoefficientFunction
            1 in the region on interest, 0 else
        """

        if self.codim != 0:
            print("Warning: Indicator is not intended for codim > 0! "
                  "Please use the IndicatorSmoothed method")
            return self.IndicatorSmoothed(lsets)

        ind_combined = CoefficientFunction(0)
        for dtt in self.as_list:
            ind = CoefficientFunction(1)

            for i, dt in enumerate(dtt):
                if dt == POS:
                    ind *= IfPos(lsets[i], 1, 0)
                elif dt == NEG:
                    ind *= IfPos(-lsets[i], 1, 0)

            ind_combined += ind
            del ind

        ind_leveled = IfPos(ind_combined, 1, 0)

        del ind_combined
        return ind_leveled

    def IndicatorSmoothed(self, lsets, eps=0.03):
        """
        Smoothed indicator function for a DomainTypeArray of codim>0.
        We assume that the level sets are approximately signed distance
        functions in the eps region of the zero set. 

        Parameters
        ----------
        lsets : tuple(ngsolve.CoefficientFunctions)
            The level set functions defining the region
        eps : float
            The distance around the subdomain which is indicated.
            (default value eps=0.01)

        Returns
        -------
        CoefficientFunction
            1 in the region of distance eps around the subdomain, 0 else
        """

        if self.codim == 0:
            print("Warning: IndicatorSmothed is not intended for codim==0! "
                  "Please use the Indicator method")
            return self.Indicator(lsets)

        def gf_abs(gf):
            return IfPos(gf, gf, -gf)

        ind_combined = CoefficientFunction(0)
        for dtt in self.as_list:
            ind = CoefficientFunction(1)
            ind_cut = CoefficientFunction(1)

            for i, dt in enumerate(dtt):
                if dt == IF:
                    ind *= IfPos(gf_abs(lsets[i]) - eps, 0, 1)
                elif dt == POS:
                    ind_cut *= IfPos(lsets[i] - eps, 1, 0)
                elif dt == NEG:
                    ind_cut *= IfPos(-lsets[i] + eps, 1, 0)

            ind_combined += ind * ind_cut
            del ind, ind_cut

        ind_leveled = IfPos(ind_combined, 1, 0)

        del ind_combined
        return ind_leveled

    def GetOuterNormals(self, lsetsp1):
        """
        For each domain region in self, we compute the outward pointing
        unit normal vector on each boundary segment of self.Boundary().

        Parameters
        ----------
        lesetsp1 : tuple(ngsolve.GridFunction)
            The set of P1 level set functions.

        Returns
        -------
        dict(dict(ngsolve.GridFunction))
            For each dtt in self, the dictionary contains a dictionary,
            where the normals can be called using the dtt of the 
            boundary segment.
        """

        if self.codim > 0:
            raise NotImplemented("GetOuterBoundary is olny available for "
                                 "codim = 0 !")

        n_dtt = {}
        bnd = self.Boundary()
        for dtt in self.as_list:
            for i, dt in enumerate(dtt):
                dtt_bnd = dtt[:i] + tuple([IF]) + dtt[i+1:]
                if dtt_bnd not in bnd:
                    continue
                if dt == POS:
                    n_dt = - 1.0 / Norm(Grad(lsetsp1[i])) * Grad(lsetsp1[i])
                else:
                    n_dt = 1.0 / Norm(Grad(lsetsp1[i])) * Grad(lsetsp1[i])
                n_dtt[dtt_bnd] = n_dt

        del bnd
        return n_dtt
