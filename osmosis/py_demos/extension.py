
from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.ngstd import *
from ngsolve.la import *

from libngsxfem_py.xfem import *
POS = 0
NEG = 1

def ExtendGF(u,fes,mesh, bfextend = None, bfmass = None, vecuold = None, setzero=False):
    outerdofs = BitArray(fes.FreeDofs())
    xfes = CastToXStdFESpace(fes)

    # no dirichlet values for the extension:
    for elid in mesh.Elements():
        for dofnr in xfes.GetDofNrs(elid):
            outerdofs.Set(dofnr)
            
    for elid in mesh.Elements():
        if xfes.XFESpace.GetDomainOfElement(elid.nr)!=POS:
            for dofnr in xfes.GetDofNrs(elid):
                outerdofs.Clear(dofnr)

    if (bfextend == None):            
        bfextend = BilinearForm(fes,"a",flags={"symmetric":True})
        bfextend.Add(CompoundBFI(bfi=BFI("laplace",coef=ConstantCF(1)), comp = 0))
        bfextend.Assemble()

    # if (bfmass == None):            
    #     bfmass = BilinearForm(fes,"a",flags={"symmetric":True})
    #     bfmass.Add(CompoundBFI(bfi=BFI("mass",coef=ConstantCF(1)), comp = 0))
    #     bfmass.Assemble()
    
    vecu = u.vec
    
    if setzero:
        for i in range(len(vecu[:])):
            if outerdofs[i]:
                if (vecuold == None):
                    vecu[i] = 0.0
                else:
                    vecu[i] = vecuold[i]
        return

    ainv = bfextend.mat.Inverse(outerdofs)
    tmp = u.vec.CreateVector()
    tmp.data = bfextend.mat * u.vec

    if (vecuold != None and bfmass != None):
        tmp.data -= bfmass.mat * vecuold 

    vecu.data -= ainv * tmp

    return

    
