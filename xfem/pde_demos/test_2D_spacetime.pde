
#
# solve the Poisson equation -Delta u = f
#
# with boundary conditions
#      u = 0  on Gamma1
#  du/dn = 1  on Gamma2

# load geometry
geometry = square.in2d

# and mesh
#mesh = square.vol.gz
#mesh = square_trigs_2x2.vol.gz
mesh = square_trigs.vol.gz
#mesh = square_quad_coarse.vol.gz

#shared = libngsxfem_test
shared = libngsxfem_xfem
#shared = libngsxfem_spacetime
#shared = libngsxfem_test

define constant heapsize = 1e7

define fespace fesh1
       -type=spacetimefes 
       -type_space=h1ho
       -order_space=1
       -order_time=1

define fespace fesx
       -type=xfespace
       -spacetime
       -levelset=(x-0.5)



numproc informxfem npix 
        -fespace=fesh1
        -xfespace=fesx
