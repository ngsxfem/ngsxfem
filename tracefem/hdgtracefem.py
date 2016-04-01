from ngsolve.comp import *
from ngsolve.fem import *
from libngsxfem_py.xfem import *
import libngsxfem_tracefem

def HDGTraceLaplaceBeltrami (coef):
    return BFI("hdgtracelaplacebeltrami", coef=coef)
