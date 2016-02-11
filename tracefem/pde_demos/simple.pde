
geometry = cube_simple.geo
#mesh = cube.vol.gz
mesh = cube_simple.vol.gz


#load xfem-library and python-bindings
shared = libngsxfem_xfem
shared = libngsxfem_tracefem


define constant heapsize = 1e9

define constant R = 1
define constant one = 1.0

####levelsets
define coefficient lset_sphere
( sqrt(x*x+y*y+z*z) - R),

define fespace fes_lset
       -type=h1ho
       -order=1

gridfunction gf_lset -fespace=fes_lset        
numproc setvalues npsv -gridfunction=gf_lset -coefficient=lset_sphere
        
define fespace fesh1
       -type=h1ho
       -order=1
       -dirichlet=[1,2,3,4]

# use an "extended" continuous finite element space
# you may change the order here
define fespace tracefes
       -type=xfespace
       -type_std=h1ho
       -ref_space=0
       -dirichlet=[1,2,3,4]

#update "extended" part of XFE space:
numproc informxfem npix
        -xfespace=tracefes
        -fespace=fesh1
        -coef_levelset=gf_lset

gridfunction u -fespace=tracefes


bilinearform a -fespace=tracefes
tracelaplacebeltrami 1.0
tracemass 1.0
#tracelaplace 0.1
# tracediv conv

bilinearform m -fespace=tracefes 
tracemass 1.0

# linearform u_zero -fespace=tracefes
# tracesource (9.0/4.0*z*z)

linearform f -fespace=tracefes
# tracesource (sin(pi*z)*(1+pi*pi*(1-z*z*z))+cos(pi*z)*4*pi*z)
#tracesource (z*z+6*z*z-2)#solution u=z*2
#tracesource (z) #solution u=z ???
tracesource (sin(pi*z)*(pi*pi*(1-z*z)+1)+cos(pi*z)*2*pi*z) #solution u=sin(pi*z)

coefficient u_sol
(sin(pi*z)),

coefficient gradu_gamma_sol
((-x*pi*z*cos(pi*z)),(-y*pi*z*cos(pi*z)),(pi*cos(pi*z)-z*pi*z*cos(pi*z))),
        
#define preconditioner c -type=local -bilinearform=a -test #-block
define preconditioner c -type=direct -bilinearform=a -inverse=pardiso -test

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6
# define preconditioner c -type=direct -bilinearform=a -inverse=pardiso -test

# numproc traceoutput npto -gridfunction=u -levelset=lset_sphere -subdivision=0 -reset -instat
# numproc traceoutput npto -gridfunction=u -levelset=lset_sphere -subdivision=0 -instat
# numproc parabolic3d np1 -bilinearforma=a -bilinearformm=m -linearform=f -visnumproc=npto -gridfunction=u -dt=0.01 -tend=3

bilinearform evalu -fespace=tracefes -nonassemble
exttrace 1.0

numproc drawflux npdf -solution=u -bilinearform=evalu -applyd -label=u

numproc visualization npviz -scalarfunction=u
    -minval=-1.5 -maxval=1.5
    -nolineartexture -deformationscale=0.25 -subdivision=0

numproc traceoutput npto -gridfunction=u -levelset=gf_lset -subdivision=0 -reset
numproc traceoutput npto -gridfunction=u -levelset=gf_lset -subdivision=0

numproc markinterface npmi -fespace=tracefes
        