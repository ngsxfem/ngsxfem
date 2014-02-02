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


define constant cuteps = 0.0

define variable t = 0.0

define coefficient lset
(y-cuteps+0.0*(t)),

numproc testxfem nptxfem 
    -levelset=lset 
    -order_space=1
    -spacetime 
    -time=t
    -timeinterval=[0,1]
    -order_time=0



