import pytest
from ngsolve import *
from xfem import *
from ngsolve.meshes import *
from netgen.geom2d import SplineGeometry
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from math import pi
from xfem.lsetcurv import *

def test_apply():
    square = SplineGeometry()
    square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=0.8))
    
    levelset = (sqrt(sqrt(x*x*x*x+y*y*y*y)) - 1.0)

    lsetp1 = GridFunction(H1(mesh))
    lsetp1.Set(levelset)
    lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}
    
    ci = CutInfo(mesh,lsetp1)
    hasneg = ci.GetElementsOfType(HASNEG)

    Vh = H1(mesh, order = 2, dirichlet=[])    
    active_dofs = GetDofsOfElements(Vh,hasneg)
    Vh = Compress(Vh,active_dofs)
    
    u,v = Vh.TnT()
    
    a = BilinearForm(Vh)
    a += SymbolicBFI(levelset_domain = lset_neg, form = u * v)
    a += SymbolicBFI(levelset_domain = lset_neg, form = grad(u)[0] * v)
    
    gfu = GridFunction(Vh)
    gfu.Set(1)

    a.Assemble()

    Au1 = gfu.vec.CreateVector()
    Au2 = gfu.vec.CreateVector()

    Au1.data = a.mat * gfu.vec
    
    a.Apply(gfu.vec,Au2)
    
    Au1.data -= Au2
    diff = Norm(Au1)
    print("diff : ",diff)
    assert diff < 1e-12

    
def test_nonlin():
    square = SplineGeometry()
    square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=0.2, quad_dominated=False))

    levelset = (sqrt(sqrt(x*x*x*x+y*y*y*y)) - 1.0)
    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lsetp1)

    Vh = H1(mesh, order=1, dirichlet=[1,2,3,4], dgjumps=True)

    ci = CutInfo(mesh, lsetp1)

    hasneg = ci.GetElementsOfType(HASNEG)  # <- "hasneg": has (also) negative level set values
    hasif = ci.GetElementsOfType(IF)
    
    ba_facets = GetFacetsWithNeighborTypes(mesh,a=hasneg,b=hasif)

    freedofs = Vh.FreeDofs()
    # freedofs &= CompoundBitArray([GetDofsOfElements(Vh,hasneg),GetDofsOfElements(Vh,haspos)])
    freedofs &= GetDofsOfElements(Vh,hasneg)

    u = Vh.TrialFunction()
    v = Vh.TestFunction()

    lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}

    a = BilinearForm(Vh, symmetric = True)
    a += SymbolicBFI(levelset_domain = lset_neg, form = (u+u**2-1)*v)
    a += SymbolicFacetPatchBFI(form = (u-u.Other())*(v-v.Other()),
                               skeleton=False,
                               definedonelements=ba_facets)

    gfu = GridFunction(Vh)

    from ngsolve.solvers import Newton

    Newton(a,gfu,freedofs=freedofs, inverse="")

    
if __name__ == "__main__":
    test_apply()
    test_nonlin()
