
geometry = cube.geo
#mesh = cube.vol.gz
mesh = cuberef2.vol.gz


#load xfem-library and python-bindings
shared = libngsxfem_xfem
shared = libngsxfem_tracefem


define constant heapsize = 1e9

define constant R = 0.666666666666666
define constant tR = 0.5
define constant tr = 0.2
define constant dtr = 0.169
define constant q = 1.21
define constant one = 1.0

####levelsets
define coefficient lset_plane
( x - R ),

define coefficient lset_sphere
( sqrt(x*x+y*y+z*z) - R),

define coefficient lset_torus
(sqrt(x*x+y*y+z*z-2*tR*sqrt(x*x+y*y)+tR*tR)-tr),

define coefficient lset_double_torus
((q*x*x+q*y*y)*(q*x*x+q*y*y)-q*x*x+q*y*y)*((q*x*x+q*y*y)*(q*x*x+q*y*y)-q*x*x+q*y*y)+z*z-dtr*dtr,

define coefficient conv
(-z*abs(z)/R,0,x*abs(z)/R)


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

#update "extended" part of XFE space:
numproc informxfem npix
        -xfespace=tracefes
        -fespace=fesh1
        -coef_levelset=lset_sphere

gridfunction u -fespace=tracefes


bilinearform a -fespace=tracefes
tracelaplacebeltrami 0.01
#tracelaplace 0.1
tracediv conv

bilinearform m -fespace=tracefes 
tracemass 1.0

linearform u_zero -fespace=tracefes
tracesource (abs(z))

linearform f -fespace=tracefes
tracesource 0


define preconditioner c -type=local -bilinearform=a -test #-block
#define preconditioner c -type=direct -bilinearform=a -inverse=pardiso -test

numproc bvp npbvp -gridfunction=u -bilinearform=m -linearform=u_zero -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6
define preconditioner c -type=direct -bilinearform=a -inverse=pardiso -test

numproc traceoutput npto -gridfunction=u -levelset=lset_sphere -subdivision=0 -reset -instat
numproc traceoutput npto -gridfunction=u -levelset=lset_sphere -subdivision=0 -instat
numproc parabolic3d np1 -bilinearforma=a -bilinearformm=m -linearform=f -visnumproc=npto -gridfunction=u -dt=0.01 -tend=3

bilinearform evalu -fespace=tracefes -nonassemble
exttrace 1.0

numproc drawflux npdf -solution=u -bilinearform=evalu -applyd -label=u

numproc visualization npviz -scalarfunction=u
    -minval=-1.5 -maxval=1.5
    -nolineartexture -deformationscale=0.25 -subdivision=0

