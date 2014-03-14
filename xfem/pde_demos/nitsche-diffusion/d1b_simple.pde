
# load geometry
geometry = d1b_simple.in2d

# and mesh
mesh = d1b_simple.vol.gz

shared = libngsxfem_xfem

define constant heapsize = 1e8

define constant zero = 0.0
define constant one = 1.0
define constant two = 2.0

define constant bneg = 1.0
define constant bpos = 1.0

define constant aneg = 1.0
define constant apos = 2.0

define constant abneg = (aneg*bneg)
define constant abpos = (apos*bpos)



define coefficient lset
#( sqrt(0.25*(x-x0)*(x-x0)+(y-y0)*(y-y0))  - R),       
(
 x
)


define coefficient sol_bnd
0,1,0,

define coefficient solneg
(2/3*(x+1)),

define coefficient solpos
(2/3+1/3*x),

define coefficient derneg
(2/3,0),

define coefficient derpos
(1/3,0),

define fespace fescomp
       -type=xh1fespace
       -order=1
       -dirichlet=[1,2]
       -empty
#       -dgjumps

numproc informxfem npix 
        -xh1fespace=fescomp
        -coef_levelset=lset

# define fespace fescomp
#        -type=h1ho
#        -order=1
#        -dirichlet=[1,2]
# #       -dgjumps

define coefficient coef_f
0,0,

define coefficient ab
((abpos)), ((abneg)), 

define gridfunction u -fespace=fescomp

define linearform f -fespace=fescomp # -print
#source coef_f -comp=1
#xsource solneg solpos

define constant small = 1e-6

define constant lambda = 5

define bilinearform a -fespace=fescomp 
       #-eliminate_internal -keep_internal -symmetric -linearform=f # -printelmat -print
#xmass one one
xlaplace aneg apos
#xnitsche_hansbo aneg apos bneg bpos lambda
#laplace ab  -comp=1


define bilinearform a_d -fespace=fescomp -nonassemble
laplace one -comp=1

numproc drawflux npdf -solution=u -bilinearform=a_d -label=grad

numproc setvalues npsv -gridfunction=u.1 -coefficient=sol_bnd -boundary

#define preconditioner c -type=local -bilinearform=a -test -block
define preconditioner c -type=direct -bilinearform=a -test
#define preconditioner c -type=bddc -bilinearform=a -test -block

##define preconditioner c -type=multigrid -bilinearform=a -test #-smoother=block

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6 # -print

numproc calccond npcc -bilinearform=a -inverse=pardiso # -gridfunction=u -printmatrix #-symmetric

numproc xdifference npxd 
        -solution=u 
        -solution_n=solneg
        -solution_p=solpos
        -derivative_n=derneg
        -derivative_p=derpos
        -levelset=lset
        -interorder=2
        -henryweight_n=(bneg)
        -henryweight_p=(bpos)


numproc visualization npviz -scalarfunction=u #-comp=0
