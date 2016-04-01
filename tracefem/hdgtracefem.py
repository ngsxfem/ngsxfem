from ngsolve.comp import *
from ngsolve.fem import *
from libngsxfem_py.xfem import *
import libngsxfem_tracefem

def HDGTraceLaplaceBeltrami (coef, param_IP_edge = 10.0, param_IP_facet = 10.0, param_normaldiffusion = 1.0):
    return BFI("hdgtracelaplacebeltrami", coef=[coef,CoefficientFunction(param_IP_edge),CoefficientFunction(param_normaldiffusion),CoefficientFunction(param_IP_facet)])
