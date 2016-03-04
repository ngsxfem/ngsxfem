from ngsolve.comp import *
from ngsolve.fem import *
from libngsxfem_py.xfem import *
import libngsxfem_tracefem

def TraceFESpace(mesh, stdfes=None, levelset=None, dgjumps=False, ref_space=0 ):
    Vh_tr = CastToXFESpace (FESpace ("xfespace", mesh=mesh, flags = {"trace" : True, "dgjumps" : dgjumps, "ref_space" : ref_space}) )
    if (stdfes!=None):
        Vh_tr.SetBaseFESpace(stdfes)
    if (levelset!=None):
        Vh_tr.SetLevelSet(levelset)
    Vh_tr.Update()
    print ("TraceFESpace-Update: done")
    return Vh_tr

def TraceMass (coef):
    return BFI("tracemass", coef=coef)

def TraceLaplaceBeltrami (coef):
    return BFI("tracelaplacebeltrami", coef=coef)

def TraceLaplace (coef):
    return BFI("tracelaplace", coef=coef)

def TraceConvection (coef):
    return BFI("tracediv", coef=coef)

def TraceSource (coef):
    return LFI("tracesource", coef=coef)

