from ngsolve.comp import *
from ngsolve.fem import *
from libngsxfem_py.xfem import *

def TwoDomainMassIntegrator (coefneg,coefpos):
    return BFI("xmass", coef=[coefneg,coefpos])

def TwoDomainSourceIntegrator (coefneg,coefpos):
    return LFI("xsource", coef=[coefneg,coefpos])



