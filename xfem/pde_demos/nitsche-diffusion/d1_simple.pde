
# load geometry
geometry = d1_simple.in2d

# and mesh
mesh = d1_simple.vol.gz

shared = libngsxfem_xfem

define constant heapsize = 1e8

# general constants
define constant zero = 0.0
define constant one = 1.0
define constant two = 2.0

# center of domain 1
define constant x0 = -0.4
define constant y0 = 0

define constant x1 = 0.4
define constant y1 = 0

define constant x2 = 0
define constant y2 = 0.4

define constant x3 = 0
define constant y3 = -0.4

# radius of domain 1
define constant R = 0.5

define coefficient lset
#( sqrt(0.25*(x-x0)*(x-x0)+(y-y0)*(y-y0))  - R),       
(
 x*x+y*y*y*y-y*y*(1-y*y)-0.1
)

numproc draw npd -coefficient=lset -label=levelset 
# henry weights
define constant bneg = 3.0
define constant bpos = 4.0

define constant aneg = 0.3
define constant apos = 0.4

define constant abneg = (aneg*bneg)
define constant abpos = (apos*bpos)


define coefficient f
(sin(x+y)),

define coefficient df
(cos(x+y)),

define coefficient mddf
(2*sin(x+y)),

define coefficient rhspos
(bpos*apos*(mddf)),

define coefficient rhsneg
(bpos*aneg*(mddf)),

define coefficient solpos
((f)),

define coefficient solposx
((df)),

define coefficient solposy
((df)),

define coefficient solneg
(bpos/bneg*(sin(x+y))),

define coefficient solnegx
(bpos/bneg*(df)),

define coefficient solnegy
(bpos/bneg*(df)),

define constant lambda = 2.0

define fespace fescomp
       -type=xh1fespace
       -order=1
       -dirichlet=[1,2]
#       -dgjumps

numproc informxfem npix 
        -xh1fespace=fescomp
        -coef_levelset=lset

define gridfunction u -fespace=fescomp

define linearform f -fespace=fescomp # -print
xsource rhsneg rhspos

define bilinearform a -fespace=fescomp #-eliminate_internal -keep_internal -symmetric -linearform=f # -printelmat -print
#xmass one one
xlaplace abneg abpos
#xnitsche_minstab_hansbo aneg apos bneg bpos
xnitsche_hansbo aneg apos bneg bpos lambda
#lo_ghostpenalty aneg apos one

numproc setvaluesx npsvx -gridfunction=u -coefficient_neg=solneg -coefficient_pos=solpos -boundary -print

#define preconditioner c -type=local -bilinearform=a -test #-block
define preconditioner c -type=direct -bilinearform=a -test
#define preconditioner c -type=bddc -bilinearform=a -test -block

##define preconditioner c -type=multigrid -bilinearform=a -test #-smoother=block

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6 # -print

define coefficient veczero
(0,0),

numproc calccond npcc -bilinearform=a -inverse=pardiso #-symmetric

numproc xdifference npxd 
        -solution=u 
        -solution_n=solneg
        -solution_p=f
        -levelset=lset
        -interorder=2
        -henryweight_n=3
        -henryweight_p=4

numproc visualization npviz -scalarfunction=u #-comp=0
