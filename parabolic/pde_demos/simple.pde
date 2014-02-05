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
mesh = square_trigs.vol.gz
#mesh = square_quad_coarse.vol.gz

#shared = libngsxfem_test
#shared = libngsxfem_common
shared = libngsxfem_parabolic

define constant heapsize = 1e7

define constant one = 1.0
define constant zero = 0.0

define constant cuteps = 0.0

define variable t = 0.0

define constant told = 0.0
define constant tnew = 0.1

define coefficient lset
( ( x - 0.2*sin(2*pi*y) > 0.7) + ( y - 0.2*sin(4*pi*x) < 0.5) + ( x - 0.2*sin(8*pi*y) < 0.3) ),
#( ( x > 0.5) * cos(2*pi*(x+y))),

define fespace fes_st 
       -type=spacetimefes 
       -type_space=h1ho
       -order_space=4
       -all_dofs_together
       -order_time=1
#       -print
#       -dirichlet=[1]

define gridfunction u_st -fespace=fes_st

define bilinearform a -fespace=fes_st # -printelmat -print
STtimeder one
STlaplace told tnew one
#STmass told tnew one
STtracepast one
#STtracefuture one
#STtracemass one one
#STtracemass zero one

define linearform f -fespace=fes_st #-print
#STsource told tnew one
STtracesource zero lset
#STtracesource one lset

numproc st_solveinstat npsi 
        -linearform=f 
        -gridfunction=u_st 
        -solver=pardiso 
        -fespace=fes_st
        -dt=5e-5
        -tend=0.05
#        -userstepping
#numproc bvp nps -bilinearform=a -linearform=f -gridfunction=u_st -solver=direct

# numproc testxfem nptxfem 
#     -levelset=lset 
#     -fespace=fes_st
#     -approx_order_space=1
#     -spacetime 
#     -time=t
#     -timeinterval=[0,1]
#     -approx_order_time=0
#     -order_time=1


define bilinearform evalu_past -fespace=fes_st 
STtracepast zero

define bilinearform evalu_future -fespace=fes_st 
STtracefuture zero

numproc drawflux npdf -bilinearform=evalu_past -solution=u_st -label=u_past
numproc drawflux npdf -bilinearform=evalu_future -solution=u_st -label=u_future

numproc visualization npvis -scalarfunction=u_future -nolineartexture -deformationscale=0.3 -subdivision=2 #-minval=0.0 -maxval=1.0
