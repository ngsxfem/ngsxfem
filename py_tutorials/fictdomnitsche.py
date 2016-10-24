from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *

from netgen.geom2d import SplineGeometry
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *

import scipy
import scipy.sparse.linalg


def SetFacetCutIndicator(mesh, gf_facet_kappa, subdivlvl=0):
    ifes = gf_facet_kappa.space
    lam = ifes.TrialFunction()
    mu = ifes.TestFunction()

    iscut = IsCut(mesh, lset_approx, subdivlvl=0)

    facetmass = BilinearForm(ifes)
    facetmass += SymbolicBFI( lam * mu, element_boundary=True)

    facetsource = LinearForm(ifes)
    facetsource += SymbolicLFI( iscut * mu, element_boundary=True)

    facetmass.Assemble()
    facetsource.Assemble()

    gf_facet_kappa.vec.data = facetmass.mat.Inverse() * facetsource.vec

    def cut(kappa):
        if (kappa > 1e-16 and kappa < 2.0-1e-16):
            return 1.0
        else:
            return 0.0
    from numpy import vectorize
    cut = vectorize(cut)
    gf_facet_kappa.vec.FV().NumPy()[:] = cut(gf_facet_kappa.vec.FV().NumPy())


from numpy import sum, sqrt as npsqrt

class IterationCounter:
    itcnt = 0
    def __init__(self, calcresfct):
        self.CalcResidual = calcresfct

    def UpdateStatus(self,v):
        # w = self.CalcResidual(v)
        self.itcnt += 1
        # print("iteration = {:3}, residual = {:0.4e}".format(self.itcnt,npsqrt(sum(w*w))))

# geometry

square = SplineGeometry()
square.AddRectangle([-.72,-.72],[.72,.72],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=2, quad_dominated=False))
mesh.Refine()
mesh.Refine()
mesh.Refine()

r = sqrt(x*x+y*y)

# sol = cos(1.5*pi*y)
# mlapsol = 2.25*pi*pi*cos(1.5*pi*y)
# levelset_ex = sqrt(y*y) - 0.33333333333333

# sol = 1.0 - 9.0 * y*y
# mlapsol = 18.0
sol=cos(2*pi*r*r)
mlapsol=8*pi*(2*pi*r*r*cos(2*pi*r*r)+sin(2*pi*r*r))
levelset_ex = r - 0.5

order = 2

lset_approx = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset_ex,lset_approx)
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2, discontinuous_qn=True)
Draw(lset_approx,mesh,"lset_approx")

Vh = H1(mesh, order=order, flags = {"dgjumps" : True})
gfu = GridFunction(Vh)

# coefficients / parameters:

n = 1.0/sqrt(InnerProduct(grad(lset_approx),grad(lset_approx))) * grad(lset_approx)
h = specialcf.mesh_size

u = Vh.TrialFunction()
v = Vh.TestFunction()

print("Vh.ndof:",Vh.ndof)

stab = 10 * order * order / h

u = Vh.TrialFunction()
v = Vh.TestFunction()


facet_cut = GridFunction(FacetFESpace(mesh,order=0))
SetFacetCutIndicator(mesh, facet_cut, subdivlvl=0)
Draw(facet_cut,mesh,"facet_cut")


lset_if  = { "levelset"       : lset_approx,
             "domain_type"    :          IF,
             "subdivlvl"      :           0}
lset_neg = { "levelset"       : lset_approx,
             "domain_type"    :         NEG,
             "subdivlvl"      :           0}

a = BilinearForm(Vh, symmetric = True, flags = { })
a += SymbolicBFI(levelset_domain = lset_neg, form = grad(u) * grad(v))

a += SymbolicBFI(levelset_domain =  lset_if,
                 form = - (grad(u) * n) * v - u * (grad(v) * n) + stab * u * v)


f = LinearForm(Vh)
f += SymbolicLFI(levelset_domain = lset_neg, form = mlapsol * v)

def dnjump(u,order):
    if order%2==0:
        return dn(u,order) - dn(u.Other(),order)
    else:
        return dn(u,order) + dn(u.Other(),order)
factors = [1.0/h,  h, h*h*h/4.0, h*h*h*h*h/9.0, h*h*h*h*h*h*h/16.0, h*h*h*h*h*h*h*h*h/25.0]
for i in range(1, order+1):
    a += SymbolicBFI( 0.2 * facet_cut * factors[i] * dnjump(u,i) * dnjump(v,i),
                      skeleton=True )

deformation = lsetmeshadap.CalcDeformation(levelset_ex)

nan = float('nan')
Draw(IfPos(-lset_approx,CoefficientFunction((0,0)),nan),mesh,"only_inner")
Draw(deformation,mesh,"deformation")
Draw(gfu,mesh,"u")

pre = Preconditioner(a,"local")

def Do(first=False):
    if not first:
        mesh.Refine()
        lset_approx.Update()
        InterpolateToP1(levelset_ex,lset_approx)
        facet_cut.Update()
        SetFacetCutIndicator(mesh, facet_cut, subdivlvl=0)
        gfu.Update()
        lsetmeshadap.CalcDeformation(levelset_ex)
    pre.Update()

    mesh.SetDeformation(deformation)
    gfu.Set(CoefficientFunction(1.0))
    a.Assemble()
    f.Assemble();



    tmp1 = f.vec.CreateVector()
    tmp2 = f.vec.CreateVector()
    def matvec(v):
        tmp1.FV().NumPy()[:] = v
        tmp2.data = a.mat * tmp1
        return tmp2.FV().NumPy()

    A = scipy.sparse.linalg.LinearOperator( (a.mat.height,a.mat.width), matvec)

    itcounter = IterationCounter(lambda x : matvec(x) - f.vec.FV().NumPy())
    gfu.vec.FV().NumPy()[:], succ = scipy.sparse.linalg.cg(
        A, f.vec.FV().NumPy(),
        callback = itcounter.UpdateStatus,
        tol=1e-8,
        maxiter=10000
    )
    print(itcounter.itcnt)
    # print(succ)
    # invmat = CGSolver(mat = a.mat, pre = pre.mat,
    #                   printrates=True, precision=1e-8, maxsteps=200)

    # gfu.vec.data = invmat * f.vec
    # last_num_its = invmat.GetSteps()
    Redraw()
    l2error = sqrt(Integrate(levelset_domain = lset_neg,
                             cf = (sol-gfu)*(sol-gfu),
                             mesh = mesh,
                             order=order))
    print("L2-Error = ", l2error)
    mesh.UnsetDeformation()
    return (l2error,itcounter.itcnt)

l2errors = []
itcnts = []

first = True

for i in range(6):
    l2error, itcnt = Do(first)
    l2errors.append(l2error)
    itcnts.append(itcnt)
    first = False
print(l2errors)
eoc = [ log(l2errors[i-1]/l2errors[i])/log(2) for i in range(1,len(l2errors))]
print(eoc)

print(itcnts)
itrate = [ log(itcnts[i-1]/itcnts[i])/log(2) for i in range(1,len(itcnts))]
print(itrate)
