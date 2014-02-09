
#
# solve the Poisson equation -Delta u = f
#
# with boundary conditions
#      u = 0  on Gamma1
#  du/dn = 1  on Gamma2

# load geometry
geometry = cube.geo

# and mesh
mesh = cube.vol.gz

#shared = libngsxfem_test
shared = libngsxfem_common
shared = libngsxfem_test

define constant heapsize = 1e7

define constant one = 1.0
define constant zero = 0.0

define constant cuteps = 0.0

define variable t = 0.0

define constant told = 0.0
define constant tnew = 0.1

define coefficient lset0
( sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))-0.233),

define coefficient lset1
( sqrt((x-0.55)*(x-0.55)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))-0.233),
#( ( x > 0.5) * cos(2*pi*(x+y))),

define fespace fes_st 
       -type=spacetimefes 
       -type_space=h1ho
       -order_space=1
#       -all_dofs_together
       -order_time=1
#       -print
#       -dirichlet=[1]

define gridfunction u_st -fespace=fes_st

define bilinearform a -fespace=fes_st # -printelmat -print
STtracepast one
STtracefuture one

define linearform f -fespace=fes_st #-print
STtracesource zero lset0
STtracesource one lset1

numproc bvp nps -bilinearform=a -linearform=f -gridfunction=u_st -solver=direct -print

numproc testxfem nptxfem 
        -levelset=u_st
        -spacetime
        -num_int_ref_space=0
        -num_int_ref_time=0

define bilinearform evalu_past -fespace=fes_st -nonassemble
STtracepast zero

define bilinearform evalu_future -fespace=fes_st -nonassemble
STtracefuture zero

numproc drawflux npdf -bilinearform=evalu_past -solution=u_st -label=u_past
numproc drawflux npdf -bilinearform=evalu_future -solution=u_st -label=u_future

numproc visualization npvis -scalarfunction=u_future -nolineartexture -deformationscale=1.0 -subdivision=2 
        -minval=-1e-6 -maxval=1e-6 -clipsolution=scalar -clipvec=[0,0,-1]
