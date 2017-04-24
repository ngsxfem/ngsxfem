#from ngsolve.fem import *
#from ngsolve.comp import *
#from ngsolve.solve import *
#from ngsolve.la import *
#from ngsolve.utils import *
#
#from xfem import *

from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *
from numpy import pi



# from ctypes import CDLL
# # on Windows replace '.so' with '.dll'
# mylngs = CDLL("libngsxfem_spacetime.so")


square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
ngmesh = square.GenerateMesh(maxh=0.05, quad_dominated=False)
mesh = Mesh (ngmesh)

fes1 = V=H1(mesh, order=1, dirichlet=[1,2,3,4])
k_t = 0
tfe = ScalarTimeFE(k_t) 

st_fes = SpaceTimeFESpace(fes1,tfe)
st_fes_ic = SpaceTimeFESpace(fes1,tfe)
st_fes.SetTime(0.5)

#visoptions.autoscale = False
#visoptions.mminval=0.0
#visoptions.mmaxval=1.0
visoptions.deformation = 1
#
#w = GridFunction(st_fes)
#
#
#Draw(w)
#
#for i in range(st_fes.ndof):
#   # print("i = {0}".format(i))
#    if i > 0:
#        w.vec[i-1] = 0
#    w.vec[i] = 1
#    Redraw()
#    #input("")
#    sleep(0.25)
    

# Fitted heat equation example (k_t = 0)
tend = 1.0
delta_t = 1/64
tnew = 0
t = Parameter(0)

trapezoidal = { "points" : [0,1], "weights" : [1/2,1/2] }

u_exact = CoefficientFunction( sin(pi*t)*sin(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y)  )
coeff_f = CoefficientFunction( pi*cos(pi*t)*sin(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y)
                               -2*pi*pi*sin(pi*t)*( cos(pi*x)*cos(pi*x)*sin(pi*y)*sin(pi*y)              
                                                   -sin(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y)
                                                   +cos(pi*y)*cos(pi*y)*sin(pi*x)*sin(pi*x)
                                                  -sin(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y))) 

u0 = GridFunction(st_fes)
u0_ic = GridFunction(st_fes_ic)
Draw(u0)
u = st_fes.TrialFunction()
v = st_fes.TestFunction()

st_fes.SetTime(1)

a = BilinearForm(st_fes,symmetric=False)
f = LinearForm(st_fes)
a.Assemble()
f.Assemble()
amat = a.mat.CreateMatrix()
fvec = f.vec.CreateVector()


while tend - tnew > delta_t/2:
              
    # clear storage
    amat.AsVector()[:] = 0
    fvec[:] = 0
    
    for ti,omega_i in zip(trapezoidal["points"],trapezoidal["weights"]):
        t.Set(tnew + delta_t*ti)
        st_fes.SetTime(ti) #for k_t = 0 this does nothing
        
        ai = BilinearForm(st_fes,symmetric=False)
        fi = LinearForm(st_fes)
        
        ai += SymbolicBFI(form = delta_t*omega_i*grad(u)*grad(v))
        fi += SymbolicLFI(form = delta_t*omega_i*coeff_f*v)
        
        if ti == 0:
            st_fes_ic.SetTime(1)
            fi += SymbolicLFI(form = u0_ic*v )
        if ti == 1:
            ai += SymbolicBFI(form = u*v )
            
    
        ai.Assemble()
        fi.Assemble()
        amat.AsVector().data += ai.mat.AsVector()
        fvec.data += fi.vec

    u0.vec.data = amat.Inverse(st_fes.FreeDofs(),"umfpack") * fvec
       
    u0_ic.vec.data = u0.vec
    
    st_fes.SetTime(1)
    tnew = tnew + delta_t
    t.Set(tnew)
    l2error = sqrt (Integrate ( (u_exact-u0)*(u_exact-u0), mesh))
           
    Redraw(blocking=True)
    
    print("t = {0}, l2error = {1}".format(tnew,l2error))
    
