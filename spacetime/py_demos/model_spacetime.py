from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *
from numpy import pi

square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
ngmesh = square.GenerateMesh(maxh=0.05, quad_dominated=False)
mesh = Mesh (ngmesh)

fes1 = V=H1(mesh, order=1, dirichlet=[1,2,3,4])
k_t = 1
tfe = ScalarTimeFE(k_t) 

st_fes = SpaceTimeFESpace(fes1,tfe)
st_fes_ic = SpaceTimeFESpace(fes1,tfe)
print("k_t = {0}".format(st_fes.k_t()))
print("Nodes of TimeFE: {0}".format(st_fes.TimeFE_nodes().NumPy()))
#visoptions.autoscale = False
#visoptions.mminval=0.0
#visoptions.mmaxval=1.0
visoptions.deformation = 1

#Fitted heat equation example
tend = 1.0
delta_t = 1/32
tnew = 0

told = Parameter(0)
tref = ReferenceTimeVariable()
t = told + delta_t*tref

trapezoidal = { "points" : [0,1], "weights" : [1/2,1/2] }
simpson = { "points" : [0,1/2,1], "weights" : [1/6,4/6,1/6] }

u_exact = CoefficientFunction( sin(pi*t)*sin(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y)  )
coeff_f = CoefficientFunction( pi*cos(pi*t)*sin(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y)
                               -2*pi*pi*sin(pi*t)*( cos(pi*x)*cos(pi*x)*sin(pi*y)*sin(pi*y)              
                                                   -sin(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y)
                                                   +cos(pi*y)*cos(pi*y)*sin(pi*x)*sin(pi*x)
                                                  -sin(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y))) 

u0 = GridFunction(st_fes)
u0_ic = GridFunction(fes1)
Draw(u0_ic)
u = st_fes.TrialFunction()
v = st_fes.TestFunction()

a = BilinearForm(st_fes,symmetric=False)
f = LinearForm(st_fes)
a.Assemble()
f.Assemble()
amat = a.mat.CreateMatrix()
fvec = f.vec.CreateVector()


# dummy lset domain to call symboliccutbfi instead of usual symbolicbfi...
levelset = (sqrt(x*x+y*y) - 1000.5)
lsetp1 = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset,lsetp1)
lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}


a = BilinearForm(st_fes,symmetric=False)
a += SymbolicBFI(levelset_domain = lset_neg, form = dt(u)*v, time_order=2)
a += SymbolicBFI(levelset_domain = lset_neg, form = delta_t*grad(u)*grad(v), time_order=2)
a += SymbolicBFI(form = fix_t(u,0)*fix_t(v,0) )
a.Assemble()

t_old = 0
u0_ic.Set(u_exact)

while tend - t_old > delta_t/2:
              
    # clear storage
    f = LinearForm(st_fes)
    f += SymbolicLFI(levelset_domain = lset_neg, form = delta_t*coeff_f*v, time_order=2)
    f += SymbolicLFI(form = u0_ic*v )
    f.Assemble()

    u0.vec.data = a.mat.Inverse(st_fes.FreeDofs(),"umfpack") * f.vec
       
    # exploiting the nodal property of the time fe:
    #u0_ic.vec[:] = u0.vec[0:fes1.ndof]
    u0_ic.vec[:] = u0.vec[fes1.ndof : 2*fes1.ndof]
    
    t_old = t_old + delta_t
    told.Set(t_old)
    
    l2error = sqrt (Integrate ( (u_exact-u0_ic)*(u_exact-u0_ic), mesh))
           
    Redraw(blocking=True)
    
    print("t = {0}, l2error = {1}".format(t_old,l2error))
    
