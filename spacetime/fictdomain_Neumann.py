# At the moment this file treats a fitted heat equation
# using the spacetimecutrule
# should evolve into an unfitted problem in the future

from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *
from numpy import pi

from xfem.lset_spacetime import *

square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
ngmesh = square.GenerateMesh(maxh=0.05, quad_dominated=False)
mesh = Mesh (ngmesh)

fes1 = V=H1(mesh, order=1, dirichlet=[1,2,3,4])
k_t = 1
k_s = 1
tfe = ScalarTimeFE(k_t) 

st_fes = SpaceTimeFESpace(fes1,tfe)
st_fes_ic = SpaceTimeFESpace(fes1,tfe)
#print("k_t = {0}".format(st_fes.k_t()))
#print("Nodes of TimeFE: {0}".format(st_fes.TimeFE_nodes().NumPy()))
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

levelset = (sqrt(x*x+y*y) - 1000.5)
lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = k_s, order_time = k_t,
                                         threshold=0.5, discontinuous_qn=True)
dfm = lset_adap_st.CalcDeformation(levelset,told,tnew,delta_t) 
lset_p1 = lset_adap_st.lset_p1

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

#lset_neg = { "levelset" : lset_p1, "domain_type" : NEG, "subdivlvl" : 0}


a = BilinearForm(st_fes,symmetric=False)
#a += SymbolicBFI(levelset_domain = lset_neg, form = dt(u)*v, time_order=1)
#a += SymbolicBFI(levelset_domain = lset_neg, form = delta_t*grad(u)*grad(v), time_order=1)
#a += SymbolicBFI(form = fix_t(u,0)*fix_t(v,0) )

a += SymbolicCutBFI(lset_p1, domain_type =  NEG, force_intorder = -1, 
                    time_order = 1, subdivlvl = 0, 
                    form = dt(u)*v )
a += SymbolicCutBFI(lset_p1, domain_type =  NEG, force_intorder = -1, 
                    time_order = 1, subdivlvl = 0, 
                    form = delta_t*grad(u)*grad(v))
a += SymbolicCutBFI(lset_p1, domain_type =  NEG, force_intorder = -1, 
                    time_order = 1, subdivlvl = 0, 
                    form = fix_t(u,0)*fix_t(v,0))


a.Assemble()

t_old = 0
u0_ic.Set(u_exact)

while tend - t_old > delta_t/2:
              
    f = LinearForm(st_fes)
    #f += SymbolicLFI(levelset_domain = lset_neg, form = delta_t*coeff_f*v, time_order=1)
    #f += SymbolicLFI(form = u0_ic*v )
    f += SymbolicCutLFI(lset = lset_p1, domain_type =  NEG, force_intorder = -1, 
                    time_order = 1, subdivlvl = 0,  form = delta_t*coeff_f*v)
    f += SymbolicCutLFI(lset = lset_p1, domain_type =  NEG, force_intorder = -1, 
                    time_order = 1, subdivlvl = 0,  form = u0_ic*v)
    f.Assemble()

    u0.vec.data = a.mat.Inverse(st_fes.FreeDofs(),"umfpack") * f.vec
       
    # exploiting the nodal property of the time fe:
    u0_ic.vec[:] = u0.vec[fes1.ndof : 2*fes1.ndof]
    
    t_old = t_old + delta_t
    told.Set(t_old)
    
    l2error = sqrt (Integrate ( (u_exact-u0_ic)*(u_exact-u0_ic), mesh))
           
    Redraw(blocking=True)
    
    print("t = {0}, l2error = {1}".format(t_old,l2error))
    

