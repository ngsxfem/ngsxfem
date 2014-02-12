
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
define constant t1 = 0.25
define constant t2 = 0.5
define constant t3 = 0.75
define constant t4 = 1.0
define constant v = 0.05

define constant cuteps = 0.0

define variable t = 0.0

define constant told = 0.0
define constant tnew = 1.0

define coefficient lset0
( sqrt((x-0.45)*(x-0.45)+(y-0.5)*(y-0.5))-1.0/sqrt(2*pi)),

define coefficient lset1
( sqrt((x-0.45-v*t1)*(x-0.45-v*t1)+(y-0.5)*(y-0.5))-1.0/sqrt(2*pi)),

define coefficient lset2
( sqrt((x-0.45-v*t2)*(x-0.45-v*t2)+(y-0.5)*(y-0.5))-1.0/sqrt(2*pi)),

define coefficient lset3
( sqrt((x-0.45-v*t3)*(x-0.45-v*t3)+(y-0.5)*(y-0.5))-1.0/sqrt(2*pi)),

define coefficient lset4
( sqrt((x-0.45-v*t4)*(x-0.45-v*t4)+(y-0.5)*(y-0.5))-1.0/sqrt(2*pi)),

# define coefficient lset0
# (x-0.25),

# define coefficient lset1
# (x-0.75),
#( ( x > 0.5) * cos(2*pi*(x+y))),

define fespace fes_st 
       -type=spacetimefes 
       -type_space=h1ho
       -order_space=1
       -all_dofs_together
       -order_time=1
#       -print
#       -dirichlet=[1]

define gridfunction u_st -fespace=fes_st

define bilinearform a -fespace=fes_st # -printelmat -print
# STtracepast one
# STtracefuture one
STtracemass t0 one
STtracemass t1 one
STtracemass t2 one
STtracemass t3 one
STtracemass t4 one

define linearform f -fespace=fes_st #-print
STtracesource t0 lset0
STtracesource t1 lset1
STtracesource t2 lset2
STtracesource t3 lset3
STtracesource t4 lset4
# STtracesource zero lset0
# STtracesource one lset1

numproc bvp nps -bilinearform=a -linearform=f -gridfunction=u_st -solver=direct -print

numproc testxfem nptxfem 
        -levelset=u_st
        -spacetime
        -num_int_ref_space=5
        -num_int_ref_time=5

define bilinearform evalu_past -fespace=fes_st -nonassemble
STtracepast zero

define bilinearform evalu_future -fespace=fes_st -nonassemble
STtracefuture zero

numproc drawflux npdf -bilinearform=evalu_past -solution=u_st -label=u_past
numproc drawflux npdf -bilinearform=evalu_future -solution=u_st -label=u_future

numproc visualization npvis -scalarfunction=u_future -nolineartexture -deformationscale=1.0 -subdivision=2 
        -minval=-1e-6 -maxval=1e-6
