from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters
from ngsolve import *
from xfem import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *

# Geometry 
square = SplineGeometry()
square.AddRectangle([-1,-1],[1,1],bc=1)

#quad version:
# ngmesh = square.GenerateMesh(maxh=1/2, quad_dominated=True)
# mesh = Mesh (ngmesh)
# mesh.Refine()
# mesh.Refine()

#trig version:
ngmesh = square.GenerateMesh(maxh=0.1, quad_dominated=False)
mesh = Mesh (ngmesh)

order = 3
# stabilization parameter for ghost-penalty
gamma_stab = [0.1,0.01,0.001,0.0001,0.00001,0.00001]
# stabilization parameter for Nitsche
lambda_nitsche  = 10 * order * order

r2 = 3/4 # outer radius
r1 = 1/4 # inner radius
rc = (r1+r2) / 2.0
rr = (r2-r1) / 2.0
r = sqrt(x*x+y*y)
levelset = IfPos(r-rc, r-rc-rr,rc-r-rr)

coeff_f = CoefficientFunction( -20*( (r1+r2)/sqrt(x*x+y*y) -4) ) 

# for monitoring the error
exact = CoefficientFunction(20*(r2-sqrt(x*x+y*y))*(sqrt(x*x+y*y)-r1))

Vh = H1(mesh, order = order, dirichlet=[],flags = {"dgjumps":True})    
gfu = GridFunction(Vh)

n_outer = specialcf.normal(mesh.dim)
h = specialcf.mesh_size   

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.1, discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)
lsetp1 = lsetmeshadap.lset_p1

lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}
lset_if  = { "levelset" : lsetp1, "domain_type" : IF , "subdivlvl" : 0}

# element, facet and dof marking w.r.t. boundary approximation with lsetp1:
ci = CutInfo(mesh,lsetp1)
hasneg = ci.GetElementsOfType(HASNEG)
active_dofs = GetDofsOfElements(Vh,hasneg)
active_dofs &= Vh.FreeDofs()

hasif = ci.GetElementsOfType(IF)
              
ba_facets = GetFacetsWithNeighborTypes(mesh,a=hasneg,b=hasif)
# facets on which ghost penalty stabilization should be applied
cf_ghost = IndicatorCF(mesh,ba_facets,facets=True)

n_levelset = 1.0/Norm(grad(lsetp1)) * grad(lsetp1)
           
a = BilinearForm(Vh,symmetric=False)
a.Assemble()
print(len(a.mat.AsVector()))
a = RestrictedBilinearForm(Vh,"test",hasneg,ba_facets)
a.Assemble()
print(len(a.mat.AsVector()))

f = LinearForm(Vh)
            
u,v = Vh.TrialFunction(), Vh.TestFunction()

# Diffusion term
diff_bfi = SymbolicBFI(lset_neg,form = grad(u)*grad(v))
diff_bfi.SetDefinedOnElements(hasneg)
a += diff_bfi
# Nitsche term
nitsche_term  = -grad(u) * n_levelset * v
nitsche_term += -grad(v) * n_levelset * u
nitsche_term += (lambda_nitsche/h) * u * v
nitsche_bfi = SymbolicBFI(lset_if,form = nitsche_term)
nitsche_bfi.SetDefinedOnElements(hasif)
a += nitsche_bfi
# rhs term:
f += SymbolicLFI(lset_neg, form=coeff_f*v)

def power(u,p):
    if p == 0: 
        return 1
    else:
        return u * power(u,p-1)

# normal derivative jumps:
def dnjump(u,order,comp = -1):
    if order%2==0:
        return dn(u,order,comp) - dn(u.Other(),order,comp)
    else:
        return dn(u,order,comp) + dn(u.Other(),order,comp)

# ghost penalty terms:
gp_term = CoefficientFunction(0.0)
for i in range(order):
    gp_term += gamma_stab[i] * power(h,2*i+1) * dnjump(u,i+1)*dnjump(v,i+1)
gp_int = SymbolicBFI(form = cf_ghost*gp_term,VOL_or_BND = VOL, skeleton=True)
gp_int.SetDefinedOnElements(ba_facets)
a += gp_int
# apply mesh adaptation    
mesh.SetDeformation(deformation)

a.Assemble()
f.Assemble()

# Solve linear system              
gfu.vec.data = a.mat.Inverse(active_dofs) * f.vec
              
#measure the error
l2error = sqrt(Integrate(lset_neg,(gfu-exact)*(gfu-exact),mesh))
print("L2 Error: {0}".format(l2error))

# unset mesh adaptation
mesh.UnsetDeformation()

#visualization:

Draw(deformation,mesh,"deformation")
Draw(levelset,mesh,"levelset")
Draw(lsetp1,mesh,"lsetp1")
Draw(gfu,mesh,"extu")
Draw(IfPos(-lsetp1,gfu,float('nan')),mesh,"u")
warped_u = CoefficientFunction((deformation[0],
                         deformation[1],
                         IfPos(-lsetp1,0.2*gfu,float('nan'))))
Draw(warped_u,mesh,"warped_u",sd=4)

from ngsolve.internal import *
visoptions.autoscale = False
visoptions.mminval=0
visoptions.mmaxval=1.25
visoptions.deformation = 1

