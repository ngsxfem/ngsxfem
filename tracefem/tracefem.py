from ngsolve.comp import *
from ngsolve.fem import *
import libngsxfem_tracefem
from xfem import *

def TraceFESpace(mesh, stdfes=None, levelset=None, dgjumps=False, ref_space=0, postpone_update = False ):
    fes = FESpace ("xfespace", mesh=mesh, flags = {"trace" : True, "dgjumps" : dgjumps, "ref_space" : ref_space})
    Vh_tr = CastToXFESpace (fes)
    if (stdfes):
        Vh_tr.SetBaseFESpace(stdfes)
    if (levelset):
        Vh_tr.SetLevelSet(levelset)
    if (postpone_update or not levelset or not stdfes):
        print ("TraceFESpace-Update: postponed")
    else:
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

def NormalLaplaceStabilization (param,normalfield):
    return BFI("normallaplacetrace", coef=[param,normalfield])

def MakeVectorToConstantFunction(fes,vec,constant=1.0):
    fesx = CastToXFESpace (fes)
    nvx = fesx.GetNVertexDofs()
    print(type(vec))
    vec[0:nvx] = 1.0
    vec[nvx:] = 0.0
