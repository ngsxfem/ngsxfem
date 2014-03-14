
# load geometry
geometry = d1_simple.in2d

# and mesh
mesh = d1_simple.vol.gz

shared = libngsxfem_xfem

define constant heapsize = 1e8

define constant zero = 0.0
define constant one = 1.0
define constant two = 2.0

define constant bneg = 1.0
define constant bpos = 1.0

define constant aneg = 1.0
define constant apos = 100.0

define constant abneg = (aneg*bneg)
define constant abpos = (apos*bpos)



define coefficient lset
#( sqrt(0.25*(x-x0)*(x-x0)+(y-y0)*(y-y0))  - R),       
(
 x*x+y*y*y*y-y*y*(1-y*y)-0.1
)

define coefficient lset_fake
(1+x),

define fespace fescomp
       -type=xh1fespace
       -order=1
       -dirichlet=[1,2]
#       -dgjumps

numproc informxfem npix 
        -xh1fespace=fescomp
        -coef_levelset=lset_fake

# define fespace fescomp
#        -type=h1ho
#        -order=1
#        -dirichlet=[1,2]
# #       -dgjumps

define coefficient coef_f
((bpos)), ((bneg)), 

define coefficient ab
((abpos)), ((abneg)), 

define gridfunction u -fespace=fescomp

define linearform f -fespace=fescomp # -print
source coef_f -comp=1

define constant small = 1e-6

define bilinearform a -fespace=fescomp 
       #-eliminate_internal -keep_internal -symmetric -linearform=f # -printelmat -print
noxlaplace aneg apos lset -comp=1
#laplace ab  -comp=1

numproc setvalues npsv -gridfunction=u.1 -coefficient=zero -boundary

define preconditioner c -type=local -bilinearform=a -test -block
#define preconditioner c -type=direct -bilinearform=a -test
#define preconditioner c -type=bddc -bilinearform=a -test -block

##define preconditioner c -type=multigrid -bilinearform=a -test #-smoother=block

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6 # -print

numproc calccond npcc -bilinearform=a -inverse=pardiso # -gridfunction=u -printmatrix #-symmetric

numproc visualization npviz -scalarfunction=u #-comp=0
