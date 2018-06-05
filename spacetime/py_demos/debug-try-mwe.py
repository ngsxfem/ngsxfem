# unfitted convection diffusion equation with Neumann b.c.
from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *
from numpy import pi

from xfem.lset_spacetime import *

# SetTestoutFile("test.out")
SetNumThreads(3)

square = SplineGeometry()
square.AddRectangle([-0.6,-0.6],[0.6,1.5],bc=1)
ngmesh = square.GenerateMesh(maxh=0.5, quad_dominated=False)
mesh = Mesh (ngmesh)

n_ref = 4
for i in range(n_ref):
    mesh.Refine()
    
h = specialcf.mesh_size

k=2
k_t = k
k_s = k
fes1  = H1(mesh,order=k_s,dirichlet=[])
fes_lset_slice = H1(mesh, order=1, dirichlet=[])
fes_dfm_slice = H1(mesh, order=k_s, dim=mesh.dim)

lset_order_time = k
time_order = 2*k_t
tfe = ScalarTimeFE(k_t) 

st_fes = SpaceTimeFESpace(fes1,tfe,flags = {"dgjumps": True, "no_low_order_space" : True})

lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = k_s, order_time = lset_order_time,
                                                    threshold=0.05, discontinuous_qn=False,periodic=False)

told = Parameter(0)
tref = ReferenceTimeVariable()


levelset= CoefficientFunction(sqrt(x*x+y*y) - 0.3)

t_old = 0.2
delta_t = 0.1
t = told + delta_t*tref
dfm = lset_adap_st.CalcDeformation(levelset,told,t_old,delta_t)
lset_p1 = lset_adap_st.lset_p1 

lset_neg = { "levelset" : lset_p1, "domain_type" : NEG, "subdivlvl" : 0}   
v = st_fes.TestFunction()  
r0 = 0.5
r1 = 0.5*pi/r0                                        
rho =  CoefficientFunction((1/(pi))*sin(2*pi*t)) 
coeff_f = CoefficientFunction(  -(pi/r0)*r1*( sin(r1*sqrt(x*x+(y-rho)*(y-rho)))*sin(r1*sqrt(x*x+(y-rho)*(y-rho))) - cos(r1*sqrt(x*x+(y-rho)*(y-rho)))*cos(r1*sqrt(x*x+(y-rho)*(y-rho))) )  + (pi/r0)*cos(r1*sqrt(x*x+(y-rho)*(y-rho)))*sin(r1*sqrt(x*x+(y-rho)*(y-rho)))*(1/sqrt(x*x+(y-rho)*(y-rho))) ) 
f = LinearForm(st_fes)
#f += SymbolicLFI(levelset_domain = lset_neg, form =  1.0*v, time_order=time_order) #works
f += SymbolicLFI(levelset_domain = lset_neg, form = coeff_f*v, time_order=time_order) #segfault

mesh.SetDeformation(dfm)
f.Assemble()
mesh.UnsetDeformation()

