from xfem import *
import libngsxfem_xstokes

def XStokesFESpace(mesh, order=1, levelset=None, dgjumps=False, empty_vel = False, dirichlet =[], ref_space=0 ):
    Vh = CastToXStokesFESpace ( FESpace ("xstokes", mesh=mesh, dirichlet=dirichlet, flags = {"order" : order, "dgjumps" : dgjumps, "ref_space" : ref_space, "dirichlet_vel" : dirichlet, "empty_vel" : empty_vel}) )
    if (levelset!=None):
        Vh.SetLevelSet(levelset)
    Vh.Update()
    print ("StokesFESpace-Update: done")
    return Vh

def TwoDomainStokesIntegrator (muneg, mupos):
    return BFI("xstokes", coef=[muneg,mupos])


def NitscheStokesBLFIntegrator (muneg, mupos, lamb=10):
    return BFI("xstokesnitsche", coef=[muneg,mupos,lamb])

def NitscheStokesLFIntegrator (gammaf):
    return LFI("xstokesnitscherhs", coef=[gammaf])

def NitscheStokesIntegrators ( muneg, mupos, lamb=10, gammaf=0 ):
    return (BFI("xstokesnitsche", coef=[muneg,mupos,lamb]), LFI("xstokesnitscherhs", coef=gammaf))


__all__ =   ["XStokesFESpace", "TwoDomainStokesIntegrator", "NitscheStokesBLFIntegrator", "NitscheStokesLFIntegrator", "NitscheStokesIntegrators"]
