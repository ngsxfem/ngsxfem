# from ngsolve.la import BaseMatrix
# from ngsolve import Projector, Norm
# from ngsolve.krylovspace import CG
from xfem import *
# from ngsolve import *
# import ngsolve
import operator
from copy import copy
from itertools import permutations, product

'''
Below we draft objects for multi level set domains
'''


class DomainTypeArray():

    def __new__(cls, dtlist):
        n = dtlist.count(ANY)
        if n == 0:
            return super(DomainTypeArray, cls).__new__(cls)
        # Else: expand ANY top POS and NEG
        dt_lists = [list(dtlist) for i in range(2**n)]
        dt_per = list(set(permutations([NEG]*n + [POS]*n, n)))

        for i, dt_list in enumerate(dt_lists):
            j = 0
            for k, dt in enumerate(dt_list):
                if dt == ANY:
                    dt_lists[i][k] = dt_per[i][j]
                    j += 1

        cdt_out = Combined_DomainType(operator.or_, DomainTypeArray(tuple(dt_lists[0])))
        for dtl in dt_lists[1:]:
            cdt_out |= DomainTypeArray(tuple(dtl))

        return cdt_out

    def __init__(self, dtlist):
        self.dtlist = dtlist
        self.codim = 0
        for m in self:
            if m == DOMAIN_TYPE.IF:
                self.codim += 1

    def __len__(self):
        return self.dtlist.__len__()

    def __contains__(self, q):
        if not type(q) in [DomainTypeArray]:
            raise Exception("need something domain-typish")
        if len(q) != len(self.dtlist):
            print("lengths don't match")
            return False
        for a, b in zip(self.dtlist, q):
            if a != b:
                return False
        return True

    def __iter__(self):
        return self.dtlist.__iter__()

    def __str__(self):
        return self.dtlist.__str__()

    def __or__(self, dta_b):
        return Combined_DomainType(operator.or_, self, dta_b)

    def __and__(self, dta_b):
        return Combined_DomainType(operator.and_, self, dta_b)

    def __invert__(self):
        return Combined_DomainType(operator.inv, self)

    def unroll(self):
        return [self.dtlist]


class Combined_DomainType():
    def __init__(self, operator, dtlist1, dtlist2=None):
        if not type(dtlist1) in [DomainTypeArray, Combined_DomainType]:
            raise Exception("need something domain-typish")
        self.dtlist1 = dtlist1
        self.codim = self.dtlist1.codim
        if dtlist2 != None:
            if not type(dtlist2) in [DomainTypeArray, Combined_DomainType]:
                raise Exception("need something domain-typish")
            self.dtlist2 = dtlist2
            if len(self.dtlist1) != len(self.dtlist2):
                raise Exception("lengths don't match")
            if self.dtlist1.codim != self.dtlist2.codim:
                raise Exception("codims don't match")
        else:
            self.dtlist2 = None
        self.operator = operator

    def __len__(self):
        return self.dtlist1.__len__()

    def __contains__(self, q):
        if self.dtlist1.codim != q.codim:
            print("codims don't match")
            return False
        a = self.dtlist1.__contains__(q)
        if self.dtlist2 != None:
            b = self.dtlist2.__contains__(q)
            return self.operator(a, b)
        else:
            if(self.operator != operator.inv):
                raise Exception("unexpected...")
            if len(self.dtlist1) != len(q):
                raise Exception("lengths don't match")
            return not a

    def __or__(self, dta_b):
        return Combined_DomainType(operator.or_, self, dta_b)

    def __and__(self, dta_b):
        return Combined_DomainType(operator.and_, self, dta_b)

    def __invert__(self):
        return Combined_DomainType(operator.inv, self)

    def __ior__(self, dta_b):
        return Combined_DomainType(operator.or_, self, dta_b)

    def __iand__(self, dta_b):
        return Combined_DomainType(operator.and_, self, dta_b)

    def __str__(self):
        return "[" + self.dtlist1.__str__() + self.operator.__str__() + self.dtlist2.__str__() + "]{codim=" + str(self.codim) + "}"

    def unroll(self):
        if type(self.dtlist1) == DomainTypeArray and type(self.dtlist2) == DomainTypeArray:
            # Combine two DomainTypeArrays
            if self.operator == operator.and_:
                dtl_out = list(self.dtlist1.dtlist)
                for i, dt in enumerate(self.dtlist1):
                    if dt != self.dtlist2.dtlist[i]:
                        if dt in [POS, NEG] and self.dtlist2.dtlist[i] in [POS, NEG]:
                            # Empty intersection
                            return []
                        elif dt == IF or self.dtlist2.dtlist[i] == IF:
                            # Increase codim, dtlist1 will be returned if non-empty intersection
                            dtl_out[i] = IF
                else:
                    return [tuple(dtl_out)]
            elif self.operator == operator.or_:
                return list(set(self.dtlist1.unroll() + self.dtlist2.unroll()))
            else:
                raise Exception("Unexpected operator (1)...")
        elif self.dtlist2 == None:
            # Only one Combined_/Domain TypeArrays to unroll
            if self.operator == operator.and_:
                return []
            elif self.operator == operator.or_:
                return self.dtlist1.unroll()
            elif self.operator == operator.inv:
                n = len(self)
                dt_per = set(permutations([NEG]*n + [POS]*n, n))
                for dt in self.dtlist1.unroll():
                    dt_per.remove(dt)
                return list(dt_per)
            else:
                raise Exception("Unexpected operator (2)...")
        elif (type(self.dtlist1) in [DomainTypeArray, Combined_DomainType] and 
              type(self.dtlist2) in [DomainTypeArray, Combined_DomainType]):
            # One Combined_TypeArray and one DomainTypeArray to unroll
            dtl1_unrolled, dtl2_unrolled = self.dtlist1.unroll(), self.dtlist2.unroll()
            if self.operator == operator.and_:
                if self.codim == 0:
                    temp = set(dtl2_unrolled) 
                    return [dtt for dtt in dtl1_unrolled if dtt in temp]
                else:
                    dtl_out = []
                    for dtt1, dtt2 in product(dtl1_unrolled, dtl2_unrolled):
                        if dtt1 == dtt2 :
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
                    return list(set(dtl_out))
            elif self.operator == operator.or_:
                return list(set(dtl1_unrolled + dtl2_unrolled))
            else:
                raise Exception("Unexpected operator (3)...")
        else:
            raise Exception("Did not expect this situation")


if __name__ == "__main__":
    print("foo")
    dtlist_a = (NEG, POS, IF, NEG)
    dta_a = DomainTypeArray(dtlist_a)
    dtlist_b = (NEG, POS, POS, NEG)
    dta_b = DomainTypeArray(dtlist_b)
    dtlist_c = (IF, IF, POS)
    dta_c = DomainTypeArray(dtlist_c)

    cdtc = ~dta_a
    print(cdtc)
    print(dta_a in cdtc)
    print(dta_b in cdtc)

    cdtc = ~dta_b
    print(cdtc)
    print(dta_a in cdtc)
    print(dta_b in cdtc)

    cdtc = ~dta_b | dta_b
    print(cdtc)
    print(dta_a in cdtc)
    print(dta_b in cdtc)

    input("Unroll tests")

    dtt_1 = (POS, NEG)
    dtt_1inv = [(NEG, POS), (POS, POS), (NEG, NEG)]
    dta_1 = DomainTypeArray(dtt_1)
    print("DomainTypeArray")
    print("Unrolled: ", dta_1.unroll())
    print("Expected: ", [dtt_1])
    input("")

    cdta_1 = Combined_DomainType(operator.or_, dta_1) 
    print("Combined_DomainType")
    print("Unrolled: ", cdta_1.unroll())
    print("Expected: ", [dtt_1])
    input("")

    cdta_2 = Combined_DomainType(operator.and_, dta_1, dta_1) 
    print("Combined_DomainType")
    print("Unrolled: ", cdta_2.unroll())
    print("Expected: ", [dtt_1])
    input("")

    cdta_3 = Combined_DomainType(operator.inv, dta_1) 
    print("Combined_DomainType: inversion")
    print("Unrolled: ", cdta_3.unroll())
    print("Expected: ", dtt_1inv)
    input("")

    dtt_2 = (ANY, POS)
    cdta_2 = DomainTypeArray(dtt_2)
    print("Combined_DomainType with ANY")
    print("Unrolled: ", cdta_2.unroll())
    print("Expected: ", [(NEG, POS), (POS, POS)])
    input("")

    cdta_3 = ~dta_1 & cdta_2
    print("Combined_DomainType inversion + and")
    print("Unrolled: ", cdta_3.unroll())
    print("Expected: ", [(NEG, POS), (POS, POS)])
    input("")

    cdta_4 = DomainTypeArray((ANY, ANY, ANY)) | DomainTypeArray((ANY, POS, NEG))
    print("Combined_DomainType ANY + (ANY + or)")
    print("Unrolled: ", cdta_4.unroll())
    print("Expected: ", list(set(permutations([NEG]*3 + [POS]*3, 3))))
    input("")

    cdta_5 = cdta_3 | DomainTypeArray((NEG, NEG))
    print("Combined_DomainType with or")
    print("Unrolled: ", cdta_5.unroll())
    print("Expected: ", [(NEG, POS), (POS, POS), (NEG, NEG)])
    input("")

    input("Test unroll for higher codim")


    cdta_6 = DomainTypeArray((POS, IF, POS)) & DomainTypeArray((IF, NEG, POS))
    print("Unrolled: ", cdta_6.unroll())
    print("Expected: ", [(IF, IF, POS)])
    input("")

    cdta_7 = DomainTypeArray((IF, POS, NEG, NEG)) | DomainTypeArray((POS, IF, NEG, NEG))
    cdta_8 = DomainTypeArray((POS, POS, IF, NEG)) | DomainTypeArray((POS, POS, NEG, IF))
    cdta_9 = cdta_7 & cdta_8

    print("Unrolled: ", cdta_9.unroll())
    print("Expected: ", [(POS, IF, IF, NEG), (IF, POS, IF, NEG),
                         (POS, IF, NEG, IF), (IF, POS, NEG, IF)])
    input("")