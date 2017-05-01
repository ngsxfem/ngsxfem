from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *
from numpy import pi


def interpol_spacetime(fes_space,lset,lset_approx,t,times):
    
    dummy_gfu = GridFunction(fes_space)
    ndof_s = len(dummy_gfu.vec)

    for i,ti in enumerate(times):
        t.Set( ti )
        dummy_gfu.Set(lset)
        lset_approx.vec[i*ndof_s : (i+1)*ndof_s] = dummy_gfu.vec[:]  
    
        
def interpolp1_spacetime(fes,fes1,lset_approx,set_approx_p1,k_t):
    
    dummy_gfu = GridFunction(fes)
    dummy_gfu_approx= GridFunction(fes1)
    ndof_s = len(dummy_gfu.vec) 
    ndof_s1 = len(dummy_gfu_approx.vec) 
    for i in range(k_t+1):
        dummy_gfu.vec[:] = lset_approx.vec[i*ndof_s : (i+1)*ndof_s]
        dummy_gfu_approx.Set(dummy_gfu)
        lset_approx_p1.vec[i*ndof_s1 : (i+1)*ndof_s1] = dummy_gfu_approx.vec[:]
           
# geometry        
square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
ngmesh = square.GenerateMesh(maxh=0.2, quad_dominated=False)
mesh = Mesh (ngmesh)

# spacetimes FEs
fes = V=H1(mesh, order=2, dirichlet=[1,2,3,4])
fes1 = V=H1(mesh, order=1, dirichlet=[1,2,3,4])
k_t = 1
tfe = ScalarTimeFE(k_t) 
st_fes = SpaceTimeFESpace(fes,tfe)
st_fes1 = SpaceTimeFESpace(fes1,tfe)

# data
t = Parameter(0)
time = 0
delta_t = 0.5
lset = CoefficientFunction( sqrt(x*x+y*y) - exp(-t)  )
lset_approx = GridFunction(st_fes)
lset_approx_p1 = GridFunction(st_fes1)

# spacetime interpolation
times = [time + delta_t * xi for xi in st_fes.TimeFE_nodes().NumPy()]
interpol_spacetime(fes,lset,lset_approx,t,times)
interpolp1_spacetime(fes,fes1,lset_approx,lset_approx_p1,k_t)

# Plotting
visoptions.deformation = 1
st_fes.SetTime(0.0)    
Draw(lset_approx)
input("")
st_fes.SetTime(1.0) 
Redraw()   
input("")
st_fes1.SetTime(0.0)    
Draw(lset_approx_p1)
input("")
st_fes1.SetTime(1.0)    
Redraw() 


