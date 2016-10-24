from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *

from netgen.geom2d import SplineGeometry
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *


def IsFacetCut(mesh,iscut, subdivlvl=0):
    ifes = FacetFESpace(mesh,order=0)
    lam = ifes.TrialFunction()
    mu = ifes.TestFunction()

    facet_kappa = GridFunction(ifes)
    facetmass = BilinearForm(ifes)
    facetmass += SymbolicBFI( lam * mu, element_boundary=True)

    facetsource = LinearForm(ifes)
    facetsource += SymbolicLFI( iscut * mu, element_boundary=True)

    facetmass.Assemble()
    facetsource.Assemble()

    facet_kappa.vec.data = facetmass.mat.Inverse() * facetsource.vec

    def cut(kappa):
        if (kappa > 1e-16 and kappa < 2.0-1e-16):
            return 1.0
        else:
            return 0.0
    from numpy import vectorize
    cut = vectorize(cut)
    facet_kappa.vec.FV().NumPy()[:] = cut(facet_kappa.vec.FV().NumPy())

    return facet_kappa


# geometry

square = SplineGeometry()
square.AddRectangle([-.72,-.72],[.72,.72],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=2, quad_dominated=False))
mesh.Refine()
mesh.Refine()
mesh.Refine()
mesh.Refine()
mesh.Refine()
mesh.Refine()

r = sqrt(x*x+y*y)
sol = 1.0 - 4.0 * r*r
mlapsol = 16.0
# sol=cos(2*pi*r*r)
# mlapsol=8*pi*(2*pi*r*r*cos(2*pi*r*r)+sin(2*pi*r*r))
levelset_ex = r - 0.5
order = 1

# levelset = GridFunction(H1(mesh,order=order+2))
# levelset.Set(levelset_ex)

lset_approx = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset_ex,lset_approx)
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2, discontinuous_qn=True)
Draw(lset_approx,mesh,"lset_approx")

Vh = H1(mesh, order=order, flags = {"dgjumps" : True})
gfu = GridFunction(Vh)

# coefficients / parameters:

n = 1.0/sqrt(InnerProduct(grad(lset_approx),grad(lset_approx))) * grad(lset_approx)
# n = 1.0/sqrt(InnerProduct(grad(levelset),grad(levelset))) * grad(levelset)
h = specialcf.mesh_size

u = Vh.TrialFunction()
v = Vh.TestFunction()

print("Vh.ndof:",Vh.ndof)

stab = 100 * order * order / h

u = Vh.TrialFunction()
v = Vh.TestFunction()

vol_cut = IsCut(mesh, lset_approx, subdivlvl=0)
facet_cut = IsFacetCut(mesh, vol_cut, subdivlvl=0)

Draw(vol_cut,mesh,"vol_cut")
Draw(facet_cut,mesh,"facet_cut")


lset_if  = { "levelset" : lset_approx, "domain_type" : IF , "subdivlvl" : 0, "force_intorder" : 2 * order}
lset_neg = { "levelset" : lset_approx, "domain_type" : NEG, "subdivlvl" : 0, "force_intorder" : 2 * order}

a = BilinearForm(Vh, symmetric = False, flags = { })
a += SymbolicBFI(levelset_domain = lset_neg, form = grad(u) * grad(v))

# a += SymbolicBFI(levelset_domain =  lset_if,
#                  form = - (grad(u) * n) * v)
a += SymbolicBFI(levelset_domain =  lset_if,
                 form = - (grad(u) * n) * v - u * (grad(v) * n) + stab * u * v)
                 # form = (grad(u)*n) * (grad(v)*n))
                 #form = - (grad(u) * n) * v + stab * u * v)
                 # form = stab * u * v)

f = LinearForm(Vh)
f += SymbolicLFI(levelset_domain = lset_neg, form = mlapsol * v)

def dnjump(u,order):
    if order%2==0:
        return dn(u,order) - dn(u.Other(),order)
    else:
        return dn(u,order) + dn(u.Other(),order)
factors = [1.0/h,  h, h*h*h, h*h*h*h*h, h*h*h*h*h*h*h]
for i in range(1, order+1):
    a += SymbolicBFI( facet_cut * factors[i] * dnjump(u,i) * dnjump(v,i), skeleton=True )

deformation = lsetmeshadap.CalcDeformation(levelset_ex)
mesh.SetDeformation(deformation)

a.Assemble()
f.Assemble();
gfu.vec.data = a.mat.Inverse(Vh.FreeDofs()) * f.vec

nan = float('nan')
# when deformation=1 and vector_function is selected, only inner/outer domain is shown:
Draw(IfPos(-lset_approx,CoefficientFunction((0,0)),nan),mesh,"only_inner")
Draw(deformation,mesh,"deformation")
Draw(gfu,mesh,"u")

print(sqrt(Integrate(levelset_domain = lset_neg,cf = (sol-gfu)*(sol-gfu),mesh = mesh,order=order)))
mesh.UnsetDeformation()
