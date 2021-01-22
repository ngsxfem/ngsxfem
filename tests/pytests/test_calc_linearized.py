import pytest
from ngsolve import *
from xfem import *
from ngsolve.meshes import *
from netgen.geom2d import SplineGeometry
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from math import pi
from xfem.lsetcurv import *

def test_calc_linearized():
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

    gfu = GridFunction(Vh)
    gfu.Set(sin(x))
    
    a1 = BilinearForm(Vh)
    a1 += SymbolicBFI(levelset_domain = lset_neg, form = u*u * v)
    a1.AssembleLinearization(gfu.vec)

    a2 = BilinearForm(Vh)
    a2 += SymbolicBFI(levelset_domain = lset_neg, form = 2*gfu*u * v)
    a2.Assemble()

    a3 = BilinearForm(Vh)
    a3 += SymbolicBFI(levelset_domain = lset_neg, form = gfu*u * v)
    
    w1 = gfu.vec.CreateVector()
    w2 = gfu.vec.CreateVector()

    w1.data = a1.mat * gfu.vec
    w2.data = a2.mat * gfu.vec
    
    w1.data -= w2
    diff = Norm(w1)
    print("diff : ",diff)
    assert diff < 1e-12

    a1.Apply(gfu.vec,w1)
    a3.Apply(gfu.vec,w2)
    
    w1.data -= w2
    diff = Norm(w1)
    print("diff : ",diff)
    assert diff < 1e-12
    
if __name__ == "__main__":
    test_calc_linearized()
    
