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
# ngmesh.Refine()
# ngmesh.Refine()
mesh = Mesh (ngmesh)

order = 2
# stabilization parameter for Nitsche
lambda_nitsche  = 10 * order * order
lambda_dg  = 10 * order * order

r2 = 3/4 # outer radius
r1 = 1/4 # inner radius
rc = (r1+r2) / 2.0
rr = (r2-r1) / 2.0
r = sqrt(x*x+y*y)
levelset = IfPos(r-rc, r-rc-rr,rc-r-rr)

coeff_f = CoefficientFunction( -20*( (r1+r2)/sqrt(x*x+y*y) -4) ) 

# for monitoring the error
exact = CoefficientFunction(20*(r2-sqrt(x*x+y*y))*(sqrt(x*x+y*y)-r1)).Compile()

Vh = L2(mesh, order = order, dirichlet=[], dgjumps = True)    
gfu = GridFunction(Vh)

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
              
ba_gp_facets = GetFacetsWithNeighborTypes(mesh,a=hasneg,b=hasif,use_and=True)
ba_fd_facets = GetFacetsWithNeighborTypes(mesh,a=hasneg,b=hasneg,use_and=True)

n_levelset = 1.0/Norm(grad(lsetp1)) * grad(lsetp1)
           
a = RestrictedBilinearForm(Vh,"a",hasneg,ba_fd_facets,check_unused=False)

f = LinearForm(Vh)
            
u,v = Vh.TrialFunction(), Vh.TestFunction()

# Diffusion term
a += SymbolicBFI(lset_neg,form = grad(u)*grad(v), definedonelements=hasneg)
a += SymbolicFacetPatchBFI(form = 0.1*1.0/h*1.0/h*(u-u.Other())*(v-v.Other()),
                           skeleton=False,
                           definedonelements=ba_fd_facets)

nF = specialcf.normal(mesh.dim)

flux_u = -0.5*(grad(u)+grad(u.Other())) * nF
flux_v = -0.5*(grad(v)+grad(v.Other())) * nF
jump_u = u-u.Other()
jump_v = v-v.Other()

#interior penalty terms:
a += SymbolicBFI(lset_neg, form = lambda_dg/h*jump_u*jump_v + flux_u * jump_v + flux_v * jump_u,skeleton=True,definedonelements=ba_fd_facets)

# Nitsche term
nitsche_term  = -grad(u) * n_levelset * v
nitsche_term += -grad(v) * n_levelset * u
nitsche_term += (lambda_nitsche/h) * u * v
a += SymbolicBFI(lset_if,form = nitsche_term, definedonelements=hasif)

# rhs term:
f += SymbolicLFI(lset_neg, form=coeff_f*v)

# apply mesh adaptation    
mesh.SetDeformation(deformation)

a.Assemble()
f.Assemble()

# Solve linear system              
gfu.vec.data = a.mat.Inverse(active_dofs,"sparsecholesky") * f.vec
              
#measure the error
l2error = sqrt(Integrate(lset_neg,(gfu-exact)*(gfu-exact),mesh))
print("L2 Error: {0}".format(l2error))

# unset mesh adaptation
mesh.UnsetDeformation()


#visualization:
import sys
if not hasattr(sys, 'argv') or len(sys.argv) == 1 or sys.argv[1] != "testmode":
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
  
  
