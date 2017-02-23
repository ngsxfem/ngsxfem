from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *

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

order = 4

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
factors = [1.0/h,  h, h*h*h/4.0, h*h*h*h*h/9.0, h*h*h*h*h*h*h/16.0, h*h*h*h*h*h*h*h*h/25.0, h*h*h*h*h*h*h*h*h*h*h/36.0]
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

    # mesh.SetDeformation(deformation)
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

    # itcounter = IterationCounter(lambda x : matvec(x) - f.vec.FV().NumPy())
    # gfu.vec.FV().NumPy()[:], succ = scipy.sparse.linalg.cg(
    #     A, f.vec.FV().NumPy(),
    #     callback = itcounter.UpdateStatus,
    #     tol=1e-8,
    #     maxiter=100000
    # )
    # print(itcounter.itcnt)

    # print(succ)
    # invmat = CGSolver(mat = a.mat, pre = pre.mat,
    #                   printrates=True, precision=1e-8, maxsteps=200)
    invmat = a.mat.Inverse()
    gfu.vec.data = invmat * f.vec
    # last_num_its = invmat.GetSteps()
    Redraw()
    l2error = sqrt(Integrate(levelset_domain = lset_neg,
                             cf = (sol-gfu)*(sol-gfu),
                             mesh = mesh,
                             order=order))
    print("L2-Error = ", l2error)
    mesh.UnsetDeformation()
    input("")
    return (l2error,1) #itcounter.itcnt)

l2errors = []
itcnts = []

first = True

for i in range(7):
    # with TaskManager():
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



# P2:
# [0.005200429999726645, 0.00044924499610759534, 5.656947260037409e-05, 6.958195622499394e-06, 6.994793922930479e-07, 9.26389405605212e-08, 1.035257049845135e-08]
# [3.533056578475361, 2.9894068085906973, 3.0232385762115195, 3.314359780232845, 2.9165908945132952, 3.1616297273023113]

# P3:
# [0.005501228805476527, 0.00013460349229927574, 5.7022939679244774e-06, 2.923598830988341e-07, 1.3864113605606791e-08, 8.124565710513994e-10]
# [5.35296616198093, 4.561029615222384, 4.285725148901028, 4.398318075287012, 4.092920873740316]

# P4:
# [0.0010816038256738505, 3.822700884929017e-05, 1.1487916617403975e-06, 3.199199532912893e-08, 6.040760886938866e-10]
# [4.822436031427234, 5.056403230026898, 5.16626239696622, 5.726836884280883]


## old (fixed eps, 2nd order FD approx. only)
# P1:
# [0.08623615939251947, 0.021613856718587357, 0.005142805920752708, 0.0012383120396572682, 0.00031075360778362615, 7.72343543252308e-05, 1.9144803269097053e-05]
# [1.9963364019580014, 2.071328910404743, 2.054180808216092, 1.9945318561938494, 2.008456525587277, 2.0122898748466502]
# [26, 49, 79, 132, 204, 375, 703]
# [-0.914270125974116, -0.6890709040618948, -0.7406133711813505, -0.6280312226130421, -0.8783214434117476, -0.9066340936892923]
# P2
# [0.005202669156860718, 0.0004486247403468851, 5.655477678473257e-05, 7.040199866401913e-06, 7.398573755038706e-07, 1.1306739247110691e-07]
# [3.5356708795969642, 2.9877883942545354, 3.0059605920071735, 3.250297296493817, 2.710064254765072]
# [275, 506, 692, 944, 1023, 1236]
# [-0.8797057662822882, -0.4516346529424145, -0.44801482172511653, -0.11594738038348747, -0.27287259815904596]
# P3
# [0.0026091988227061207, 0.00013316646178441238, 5.790494126948221e-06, 3.8491593946945305e-07, 4.452912171052677e-08]
# [4.292304194778195, 4.523400507924929, 3.911071144759691, 3.1117223532119014]
# [1221, 2170, 3629, 4533, 5396]
# [-0.8296318423444384, -0.7418770141474926, -0.32089410440665106, -0.2514241870622823]
