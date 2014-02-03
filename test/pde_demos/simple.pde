#
# solve the Poisson equation -Delta u = f
#
# with boundary conditions
#      u = 0  on Gamma1
#  du/dn = 1  on Gamma2

# load geometry
geometry = square.in2d

# and mesh
mesh = square.vol.gz

shared = libngsxfem_test
shared = libngsxfem_common


define constant cuteps = 0.0

define variable t = 0.0

define coefficient lset
(y-cuteps+0.0*(t)),

define fespace fes_st 
       -type=spacetimefes 
       -type_space=l2ho
       -order_space=1 
       -order_time=1 
       -print
       -dirichlet=[0,1]

define gridfunction u_st -fespace=fes_st

numproc testxfem nptxfem 
    -levelset=lset 
    -fespace=fes_st
    -approx_order_space=1
    -spacetime 
    -time=t
    -timeinterval=[0,1]
    -approx_order_time=0
    -order_time=1



