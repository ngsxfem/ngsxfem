from ngsolve.comp import *
from ngsolve.fem import *
from libngsxfem_py.xfem import *
import libngsxfem_tracefem

def TraceFESpace(mesh, stdfes=None, levelset=None, dgjumps=False, ref_space=0, postpone_update = False ):
    fes = FESpace ("xfespace", mesh=mesh, flags = {"trace" : True, "dgjumps" : dgjumps, "ref_space" : ref_space})
    Vh_tr = CastToXFESpace (fes)
    if (stdfes!=None):
        Vh_tr.SetBaseFESpace(stdfes)
    if (levelset!=None):
        Vh_tr.SetLevelSet(levelset)
    if (postpone_update or levelset==None or stdfes==None):
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
    vec[:][0:nvx] = 1.0
    vec[:][nvx:] = 0.0
