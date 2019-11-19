# from ngsolve.la import BaseMatrix
# from ngsolve import Projector, Norm
# from ngsolve.krylovspace import CG
from xfem import *
# from ngsolve import *
# import ngsolve
import operator
from copy import copy

'''
Below we draft objects for multi level set domains
'''

class DomainTypeArray():
  def __init__(self,dtlist):
    self.dtlist = dtlist      
    self.codim = 0
    for m in self:
      if m == DOMAIN_TYPE.IF:
        self.codim += 1

  def __len__(self):
    return self.dtlist.__len__()
    
  def __contains__(self,q):
    if not type(q) in [DomainTypeArray]:
       raise Exception("need something domain-typish")
    if len(q) != len(self.dtlist):
        print("lengths don't match")
        return False
    for a,b in zip(self.dtlist,q):
        if a!=b:
            return False
    return True

  def __iter__(self):
    return self.dtlist.__iter__()
  def __str__(self):
    return self.dtlist.__str__()

  def __or__(self,dta_b):
    return Combined_DomainType(operator.or_,self,dta_b)
  def __and__(self,dta_b):
    return Combined_DomainType(operator.and_,self,dta_b)
  def __invert__(self):
    return Combined_DomainType(operator.inv,self)


    
class Combined_DomainType():
  def __init__(self,operator,dtlist1,dtlist2=None):
    if not type(dtlist1) in [DomainTypeArray,Combined_DomainType]:
       raise Exception("need something domain-typish")
    self.dtlist1 = dtlist1
    self.codim = self.dtlist1.codim
    if dtlist2 != None:
      if not type(dtlist2) in [DomainTypeArray,Combined_DomainType]:
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

  def __contains__(self,q):
    if self.dtlist1.codim != q.codim:
      print("codims don't match")
      return False
    a = self.dtlist1.__contains__(q)
    if self.dtlist2 != None:
      b = self.dtlist2.__contains__(q)
      return self.operator(a,b)
    else:
      if(self.operator != operator.inv):
        raise Exception("unexpected...")
      if len(self.dtlist1) != len(q):
        raise Exception("lengths don't match")
      return not a

  def __or__(self,dta_b):
    return Combined_DomainType(operator.or_,self,dta_b)
  def __and__(self,dta_b):
    return Combined_DomainType(operator.and_,self,dta_b)
  def __invert__(self):
    return Combined_DomainType(operator.inv,self)

  def __str__(self):
    return "[" + self.dtlist1.__str__()+self.operator.__str__()+self.dtlist2.__str__() + "]{codim=" + str(self.codim) + "}"
    
if __name__ == "__main__":
  print("foo")
  dtlist_a = (NEG,POS,IF,NEG)
  dta_a = DomainTypeArray(dtlist_a)
  dtlist_b = (NEG,POS,POS,NEG)
  dta_b = DomainTypeArray(dtlist_b)
  dtlist_c = (IF,IF,POS)
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
