
# load geometry
geometry = d1a_simple.in2d

# and mesh
mesh = d1a_simple6.vol.gz

shared = libngsxfem_xfem
shared = libngsxfem_parabolic

define constant heapsize = 1e9

define constant zero = 0.0
define constant one = 1.0
define constant two = 2.0

define constant bneg = 2.0
define constant bpos = 1.0

define constant aneg = 1.0
define constant apos = 5.0

define constant abneg = (aneg*bneg)
define constant abpos = (apos*bpos)

define constant aobneg = (aneg/bneg)
define constant aobpos = (apos/bpos)

define constant R = 0.5+1e-10


define coefficient lset
( sqrt(x*x+y*y) - R),       

define coefficient solneg
(apos*cos((x*x+y*y)/(R*R)*pi/2)+bpos),

define coefficient solpos
(aneg*cos((x*x+y*y)/(R*R)*pi/2)+bneg),

define coefficient rhsneg
(aneg*apos*bneg*pi/(R*R)*(2*sin((x*x+y*y)/(R*R)*pi/2)+pi/(R*R)*(x*x+y*y)*cos((x*x+y*y)/(R*R)*pi/2))),

define coefficient rhspos
(aneg*apos*bpos*pi/(R*R)*(2*sin((x*x+y*y)/(R*R)*pi/2)+pi/(R*R)*(x*x+y*y)*cos((x*x+y*y)/(R*R)*pi/2))),

define fespace fescomp
       -type=xh1fespace
       -order=2
       -dirichlet=[1,2,3,4]
       -ref_space=8
#       -empty
#       -dgjumps

numproc informxfem npix 
        -xh1fespace=fescomp
        -coef_levelset=lset

# define fespace fescomp
#        -type=h1ho
#        -order=1
#        -dirichlet=[1,2]
# #       -dgjumps

# define coefficient coef_f
# 0,0,

# define coefficient ab
# ((abpos)), ((abneg)), 

define gridfunction u -fespace=fescomp

define linearform f -fespace=fescomp # -print
xsource rhsneg rhspos
#source coef_f -comp=1
#xsource solneg solpos

define constant small = 1e-6

define constant lambda = 2

define bilinearform a -fespace=fescomp 
#-eliminate_internal -keep_internal -symmetric -linearform=f # -printelmat -print
#xmass one one
xlaplace abneg abpos
xnitsche_hansbo aneg apos bneg bpos lambda
#laplace ab  -comp=1


define bilinearform a_d -fespace=fescomp -nonassemble
laplace one -comp=1

numproc drawflux npdf -solution=u -bilinearform=a_d -label=grad

# numproc setvalues npsv -gridfunction=u.1 -coefficient=solneg -boundary

numproc setvaluesx npsvx -gridfunction=u -coefficient_neg=solneg -coefficient_pos=solpos -boundary # -print

#define preconditioner c -type=local -bilinearform=a #-test #-block
define preconditioner c -type=direct -bilinearform=a -inverse=pardiso #-test
#define preconditioner c -type=bddc -bilinearform=a -test -block

#define preconditioner c -type=spacetime -bilinearform=a #-test #-smoother=block

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=100000 -prec=1e-6 # -print

# numproc calccond npcc -bilinearform=a -inverse=pardiso # -gridfunction=u -printmatrix #-symmetric

numproc xdifference npxd 
        -solution=u 
        -solution_n=solneg
        -solution_p=solpos
        -derivative_n=derneg
        -derivative_p=derpos
        -levelset=lset
        -interorder=5
        -henryweight_n=2
        -henryweight_p=1
        -diffusion_n=1.0
        -diffusion_p=5.0


numproc visualization npviz -scalarfunction=u #-comp=0
