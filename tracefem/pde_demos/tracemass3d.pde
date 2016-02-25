
# geometry = square.in2d
# mesh = square.vol.gz
geometry = cube.geo
mesh = cube.vol.gz

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

# interface description as zero-level

define coefficient lset_plane
( x - R ),

define coefficient lset_sphere
( sqrt(x*x+y*y+z*z) - R),

define coefficient lset_torus
(sqrt(x*x+y*y+z*z-2*tR*sqrt(x*x+y*y)+tR*tR)-tr),

define coefficient lset_double_torus
((q*x*x+q*y*y)*(q*x*x+q*y*y)-q*x*x+q*y*y)*((q*x*x+q*y*y)*(q*x*x+q*y*y)-q*x*x+q*y*y)+z*z-dtr*dtr,

define fespace fesh1
       -type=h1ho
       -order=1

# use an "extended" continuous finite element space
# you may change the order here
define fespace tracefes
       -type=xfespace
       -type_std=h1ho
#       -ref_space=1
#       -empty

#update "extended" part of XFE space:
numproc informxfem npix
        -xfespace=tracefes
        -fespace=fesh1
        # -coef_levelset=lset_plane
        # -coef_levelset=lset_sphere
        -coef_levelset=lset_torus
        # -coef_levelset=lset_double_torus

gridfunction u -fespace=tracefes

bilinearform a -fespace=tracefes -symmetric
tracemass 1.0

linearform f -fespace=tracefes
tracesource (1+sin(pi*2*x*y*z))

define preconditioner c -type=local -bilinearform=a #-test #-block
#define preconditioner c -type=direct -bilinearform=a -inverse=pardiso #-test

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6

bilinearform evalu -fespace=tracefes -nonassemble
exttrace 1.0

numproc drawflux npdf -solution=u -bilinearform=evalu -applyd -label=u

#numproc draw npdf2 -coefficient=lset_plane -label=levelset_plane
#numproc draw npdf3 -coefficient=lset_sphere -label=levelset_sphere
#numproc draw npdf4 -coefficient=lset_torus -label=levelset_torus
#numproc draw npdf5 -coefficient=lset_double_torus -label=levelset_doubllge_torus

numproc visualization npviz 
    -clipvec=[-1,-1,-1]
    -centerpoint=[0,0,0]
    -clipsolution=scalar
    -scalarfunction=u
    -minval=0.5 -maxval=1.5
    -nolineartexture -deformationscale=0.25 -subdivision=0

numproc markinterface npmi -fespace=tracefes

numproc traceoutput npto -gridfunction=u -levelset=lset_torus -subdivision=0

#TODO visualize:
# * use xfem visualizer (done in the extension sense...)
# or
# * write your own visualization routine
