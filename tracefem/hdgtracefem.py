from ngsolve.comp import *
from ngsolve.fem import *
from libngsxfem_py.xfem import *
import libngsxfem_tracefem

def HDGTraceLaplaceBeltrami (coef, param_lambda = 10.0 ):
    return BFI("hdgtracelaplacebeltrami", coef=[coef,CoefficientFunction(param_lambda)])
