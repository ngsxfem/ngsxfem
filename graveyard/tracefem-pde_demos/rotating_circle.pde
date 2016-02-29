
# load geometry
geometry = square2.in2d
# and mesh
mesh = square2.vol.gz

#load xfem-library and python-bindings
shared = libngsxfem_xfem
shared = libngsxfem_tracefem

define constant heapsize = 1e9

define constant R = 1#0.333333333333333333
define constant one = 1.0

# interface description as zero-level
define coefficient lset
#( x - R ),
( sqrt(x*x+y*y) - R),

define coefficient conv
#(-100*y*y*x/(x*x+y*y),-100*y*y*y/(x*x+y*y))
#(2*y*y*y*x/((x*x+y*y)*(x*x+y*y)*sqrt(x*x+y*y)),-2*y*y*x*x/((x*x+y*y)*(x*x+y*y)*sqrt(x*x+y*y))
(-y,x)
define fespace fesh1
       -type=h1ho
       -order=1
       -dirichlet=[1,2,3,4]

# use an "extended" continuous finite element space
# you may change the order here
define fespace tracefes
       -type=xfespace
       -type_std=h1ho
       -ref_space=1
        -dirichlet=[1,2,3,4]
#       -empty

#update "extended" part of XFE space:
numproc informxfem npix
        -xfespace=tracefes
        -fespace=fesh1
        -coef_levelset=lset

gridfunction u -fespace=tracefes


bilinearform a -fespace=tracefes
#tracelaplacebeltrami 0.1
tracelaplace 0.1
tracediv conv

bilinearform m -fespace=tracefes 
tracemass 1.0

linearform u_zero -fespace=tracefes
tracesource (x)


linearform f -fespace=tracefes
#tracesource sin(pi*y/1.5)
tracesource 0




define preconditioner c -type=local -bilinearform=a -test #-block
#define preconditioner c -type=direct -bilinearform=a -inverse=pardiso -test

numproc bvp npbvp -gridfunction=u -bilinearform=m -linearform=u_zero -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6
define preconditioner c -type=direct -bilinearform=a -inverse=pardiso -test


numproc parabolic np1 -bilinearforma=a -bilinearformm=m -gridfunction=u -linearform=f  -dt=0.0001 -tend=10

bilinearform evalu -fespace=tracefes -nonassemble
exttrace 1.0

numproc drawflux npdf -solution=u -bilinearform=evalu -applyd -label=u

#numproc draw npdf2 -coefficient=lset -label=levelset

numproc visualization npviz -scalarfunction=u
    -minval=-1.5 -maxval=1.5
    -nolineartexture -deformationscale=0.25 -subdivision=0

#numproc markinterface npmi -fespace=tracefes


#TODO visualize:
# * use xfem visualizer (done in the extension sense...)
# or
# * write your own visualization routine
