"""
Convenience layer module for integration using multiple level sets.
"""
from ngsolve import Norm, Grad, GridFunction, CoefficientFunction, IfPos
from xfem import IF, POS, NEG, ANY, Integrate
from itertools import chain, product, repeat
from collections import Counter


class DomainTypeArray():
    """
    Container class for multiple level sets.

    Parameters
    ----------
    dtlist : {list(tuples({ENUM,ANY})), tuple({ENUM,ANY})}
        The tuple(s) describing the region of interest. The regions
        described by the individual tuples must be of the same
        co-dimension and the same length. Each ANY is expanded
        internally into a POS and a NEG region.
    lsets : tuple(ngsolve.GridFunction)
        Optional argument (default=None): A tuple of discrete level set
        functions, with respect to which the domain in defined. If this
        given, then the DomainTypeArray is compressed into regions with
        positive measure on initialisation.
    persistent_compress : bool
        Optional argument (default=False): Sets whether the result from
        operations on this DomainTypeArray instance should automatically
        be compressed with respect to lsets.

    Attributes
    ----------
    as_list : list(tuples)
        Expanded list containing tuples describing the region of
        interest.
    codim : int
        Co-dimension of the contained region of interest.
    lsets : {None, tuple(ngsolve.GridFunction)}
        The set of level set functions, with respect to which the
        instance is compressed if persistent_compress=True.
    persistent_compress : bool
        Indicates, whether the instance to be compressed persistently
        after operations. the default value is False.
    do_safety_checks : bool
        Option to turn of safety checks on the input list. This is
        useful if a copy of an existing DomainTypeArray is made and it
        is known, that the input is permissible. Default is value is
        True.

    Methods
    -------
    Boundary(element_marking=False):
        A DomainTypeArray describing the boundary of the current
        DomainTypeArray region.
    Compress(lsets, persistent=False):
        Remove domain regions which have have zero measure with respect
        to the given level set functions.
    Indicator(lsets):
        Returns an indicator CoefficientFunction of the instances'
        region (codim=0 only).
    IndicatorSmoothed(lsets, eps):
        Returns an indicator CoefficientFunction of an eps-region 
        around the current instances' region (codim>0 only).
    GetOuterNormals(lsetsp1):
        Returns a dictionary, containing the outward pointing unit
        normal vectors on each section of self.Boundary(). The keys are
        the tuples of the boundary segments, each normal vector is
        defined on.
    """

    def __init__(self, dtlist, lsets=None, persistent_compress=False,
                 do_safety_checks=True):
        if type(dtlist) == tuple:
            dtlist = [dtlist]
        if type(dtlist) != list or type(dtlist[0]) != tuple:
            raise TypeError("Invalid input: dtlist must be a tuple or a list of tuples")

        dtt_len = len(dtlist[0])

        self.codim = sum([1 for dt in dtlist[0] if dt == IF])
        if self.codim > 3:
            print("Warning: Codim > 3 !!!")
        if do_safety_checks:
            for dtt in dtlist:
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
            if lsets:
                if len(lsets) != dtt_len:
                    raise Exception("The number of level sets does not match the "
                                    " length of the arrays domain tuples!")
                if type(lsets[0]) != GridFunction:
                    raise Exception("The level set functions need to be "
                                    "ngsolve.GridFunctions!")

        # Expansion of ANY
        if (ANY in chain.from_iterable(dtlist)):
            dtlist_expand = []
            for dtt in dtlist:
                if ANY in dtt:
                    locations_expand = [i for i in range(len(dtt)) if dtt[i] == ANY]
                    n = len(locations_expand)
                    dt_lists = [list(dtt) for i in range(2**n)]
                    dt_per = list(product(*list(repeat((POS, NEG), n))))
                    for i, dtl in enumerate(dt_lists):
                        for j, k in enumerate(locations_expand):
                            dt_lists[i][k] = dt_per[i][j]
                        dtlist_expand.append(tuple(dt_lists[i]))
                else:
                    dtlist_expand.append(dtt)
            self.as_list = list(set(dtlist_expand))
        else:
            self.as_list = list(set(dtlist))

        self.lsets = lsets
        self.persistent_compress = persistent_compress

        if self.lsets:
            self.Compress(self.lsets, self.persistent_compress)
        
        if not self.persistent_compress:
            self.lsets = None

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
        if not self.persistent_compress and not dta_b.persistent_compress:
            lsets_out = None
            pers_out = False
        elif ((self.persistent_compress and  not dta_b.persistent_compress) or 
              (not self.persistent_compress and dta_b.persistent_compress)):
            print("WARNING: Combined domains with non matching persistent "
                  "compression status: \n The returned domain will not be "
                  "compressed.")
            lsets_out = None
            pers_out = False
        else:
            for lset1, lset2 in zip(self.lsets, dta_b.lsets):
                if lset1 != lset2:
                    raise Exception("Operator only possible for domains with the same level sets")
            lsets_out = self.lsets
            pers_out = self.persistent_compress

        return DomainTypeArray(self.as_list + dta_b.as_list, lsets_out, pers_out)

    def __and__(self, dta_b):
        if not self.persistent_compress and not dta_b.persistent_compress:
            lsets_out = None
            pers_out = False
        elif ((self.persistent_compress and  not dta_b.persistent_compress) or 
              (not self.persistent_compress and dta_b.persistent_compress)):
            print("WARNING: Combined domains with non matching persistent "
                  "compression status: \n The returned domain will not be "
                  "compressed.")
            lsets_out = None
            pers_out = False
        else:
            for lset1, lset2 in zip(self.lsets, dta_b.lsets):
                if lset1 != lset2:
                    raise Exception("Operator only possible for domains with the same level sets")
            lsets_out = self.lsets
            pers_out = self.persistent_compress

        dtl_out = [dtt for dtt in self.as_list if dtt in dta_b.as_list]
        
        if len(dtl_out) == 0:
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
        
        return DomainTypeArray(dtl_out, lsets_out, pers_out)

    def __invert__(self):
        n = len(self.as_list[0])
        if self.codim == 0:
            dt_per = set(product(*list(repeat((POS, NEG), n))))
        else:
            dt_per = set(product(*list(repeat((POS, NEG, IF), n))))
            dt_per = set(filter(lambda dtt: dtt.count(IF) == self.codim, dt_per))
        dt_per = list(dt_per.difference(set(self.as_list)))
        return DomainTypeArray(dt_per, self.lsets, self.persistent_compress)

    def __ior__(self, dta_b):
        return self.__or__(dta_b)

    def __iand__(self, dta_b):
        return self.__and__(dta_b)

    def Compress(self, lsets, persistent=False):
        """
        For a given set of level sets, this function removes any regions
        contained in the instance, which do not have a positive measure.
        The level sets remain attached to the instance, such that 
        self.Boundary() called after self.Compress(lsets) automatically
        compresses the output. 
        
        Parameters
        ----------
        lsets : tuple(ngsolve.GridFunction)
            The set of level sets with respect to which it is computed,
            whether regions have positive measure.
        persistent : bool
            Sets compression to continue on operations such as 
            Boundary(). Default vale = False 

        Returns
        -------
        None
        """

        if len(lsets) != len(self.as_list[0]):
            raise Exception("The number of level sets does not match the "
                            " length of the arrays domain tuples!")
        if type(lsets[0]) != GridFunction:
            raise Exception("The level set functions need to be "
                            "ngsolve.GridFunctions!")

        mesh = lsets[0].space.mesh
        dtl_out = []
        for dtt in self.as_list:
            weight = Integrate({"levelset": lsets, "domain_type": dtt},
                               cf=1, mesh=mesh, order=0)
            if abs(weight) > 1e-12:
                dtl_out.append(dtt)
        self.as_list = dtl_out

        if persistent:
            self.persistent_compress = persistent
            self.lsets = lsets

        return None


    def Boundary(self):
        """
        Compute the boundary of the instances' contained region. If
        persistent_compress=True has been set for the instance, then
        the result will compressed into boundary regions with positive
        measure. 

        Returns
        -------
        DomainTypeArray
            Container of the resulting boundary region.
        """
        if self.codim >= 3:
            raise Exception("Boundary does not make sense for codim >=3")

        dtl_out = []
        for dtt in self.as_list:
            for i, dt in enumerate(dtt):
                if dt != IF:
                    if (dtt[:i] + (POS,) + dtt[i+1:] in self and
                        dtt[:i] + (NEG,) + dtt[i+1:] in self):
                        continue
                    else:
                        is_relevant = True
                        dtt_out = dtt[:i] + (IF,) + dtt[i+1:]
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

        if self.persistent_compress:
            dta_out = DomainTypeArray(dtl_out, self.lsets, 
                                      self.persistent_compress)
        else:
            dta_out = DomainTypeArray(dtl_out)
    
        return dta_out

    def Indicator(self, lsets):
        """
        Indicator CoefficientFunction for a DomainTypeArray of codim=0.

        Parameters
        ----------
        lsets : tuple(ngsolve.CoefficientFunctions)
            The set of level set functions with respect to which the
            instances' region is defined.

        Returns
        -------
        CoefficientFunction
            1 in the region on interest, 0 else.
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
        Indicator CoefficientFunction, indicating the eps-region around
        the instances' region (codim>0). We assume that the level sets
        are approximately signed distance functions in the eps region of
        the zero set. 

        Parameters
        ----------
        lsets : tuple(ngsolve.CoefficientFunctions)
            The set of level set functions defining the region.
        eps : float
            The distance around the sub-domain which is indicated
            (default eps=0.01).

        Returns
        -------
        CoefficientFunction
            1 in the region of distance eps around the sub-domain, 0 
            else.
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
        For each region in self.Boundary(), we compute the outward 
        pointing unit normal vector.

        Parameters
        ----------
        lesetsp1 : tuple(ngsolve.GridFunction)
            The set of P1 level set functions.

        Returns
        -------
        dict(ngsolve.GridFunction)
            Each normal is accessible with the tuple defining the
            boundary region it is defined on as the key.
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


def TensorUnion(*args):
    """
    Construct the union of an arbitrary number of regions, described by
    different DomainTypeArrays, each dependent on a different set of 
    level sets. This increases the tuple-dimension of the resulting 
    domain regions.

    Parameters
    ----------
    *args : DomainTypeArray
        Variable list of DomainTypeArray objects of the same
        co-dimension.

    Returns
    -------
    DomainTypeArray
        Container of the resulting region.
    """

    n_dtas = []
    codim = args[0].codim

    for dta in args:
        if type(dta) != DomainTypeArray:
            raise TypeError("TensorUnion only possible for DomainTypeArrays")
        if dta.codim != codim:
            raise Exception("Cannot form TensorUnion for arrays of different codimension")
        n_dtas.append(len(dta.as_list[0]))

    i = 0
    dtl_out = []
    dtt_candiate = [ANY for j in range(sum(n_dtas))]
    for j, dta in enumerate(args):
        for dtt in dta.as_list:
            dtt_out = dtt_candiate.copy()
            dtt_out[i:i+n_dtas[j]] = dtt   
            dtl_out.append(tuple(dtt_out))
        i += n_dtas[j]

    return DomainTypeArray(dtl_out)


def TensorIntersection(*args):
    """
    Construct the intersection of an arbitrary number of regions, 
    described by different DomainTypeArrays, each dependent on a
    different set of level sets. This increases the tuple-dimension of
    the resulting domain regions.

    Parameters
    ----------
    *args : DomainTypeArray
        Variable list of DomainTypeArray objects.

    Returns
    -------
    DomainTypeArray
        Container of the resulting region.
    """
    
    for dta in args:
        if type(dta) != DomainTypeArray:
            raise TypeError("TensorUnion only possible for DomainTypeArrays")

    dtl_out = []
    for dtt in product(*[dta.as_list for dta in args]):
        dtl_out.append(tuple(chain(*dtt)))
    
    return DomainTypeArray(dtl_out)
