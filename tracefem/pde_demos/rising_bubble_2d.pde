
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

define constant v0=20
define constant mu_in=10
define constant mu_out=1

define constant conv_coeff=1/2*v0*mu_out/(mu_in+mu_out)



# interface description as zero-level
define coefficient lset
#( x - R ),
( sqrt(x*x+y*y) - R),

define coefficient conv
#(-100*y*y*x/(x*x+y*y),-100*y*y*y/(x*x+y*y))
(conv_coeff/(R*R)*y*x,-conv_coeff/(R*R)*x*x)

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
#tracelaplacebeltrami 0.001
tracelaplace 0.001
tracediv conv

bilinearform m -fespace=tracefes 
tracemass 1.0

linearform u_zero -fespace=tracefes
#tracesource exp(-10*(x-1)*(x-1))
#tracesource (2*sin(10*atan(y/x))+2)*(y+1)
#tracesource 1

linearform f -fespace=tracefes
#tracesource sin(pi*y/1.5)
#tracesource 5*exp(-10000*(x+1)*(x+1)*(x+1)*(x+1))
#tracesource 0
tracesource 10*exp(-100*((y-1)*(y-1)))#+(x-sqrt(R*R-0.8*0.8))*(x-sqrt(R*R-0.8*0.8))))
#tracesource cos(pi/2*x)*y/abs(y)#(y+abs(y))/y


define preconditioner c -type=local -bilinearform=a -test #-block
#define preconditioner c -type=direct -bilinearform=a -inverse=pardiso -test

numproc bvp npbvp -gridfunction=u -bilinearform=m -linearform=u_zero -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6
define preconditioner c -type=direct -bilinearform=a -inverse=pardiso -test


numproc parabolic np1 -bilinearforma=a -bilinearformm=m -gridfunction=u -linearform=f  -dt=0.0001 -tend=10 -periodic_rhs=3

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
