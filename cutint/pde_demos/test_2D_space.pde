
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
shared = libngsxfem_common
shared = libngsxfem_test

define constant heapsize = 1e7

define constant one = 1.0
define constant zero = 0.0
define constant t0 = 0.0
define constant v = 0.0
define constant s = 0.0

define constant cuteps = 0.0

define variable t = 0.0

define coefficient lset
( sqrt((x-0.45)*(x-0.45)+(y-0.5)*(y-0.5))-(1.0-t0*s)/sqrt(2*pi)),

define fespace fes 
       -type=l2ho
       -order=3
       -all_dofs_together

define gridfunction u -fespace=fes

define bilinearform a -fespace=fes # -printelmat -print
mass one 

define linearform f -fespace=fes #-print
source lset

numproc bvp nps -bilinearform=a -linearform=f -gridfunction=u -solver=direct -print

numproc testxfem nptxfem 
        -levelset=u
        -num_int_ref_space=3
        -int_order_space=2
        -bound
        

numproc visualization npvis -scalarfunction=u -nolineartexture -deformationscale=1.0 -subdivision=2 
        -minval=-1e-6 -maxval=1e-6
