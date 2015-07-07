
# load geometry
geometry = d2_xnitsche.in2d

# and mesh
mesh = d2_xnitsche.vol.gz


shared = libngsxfem_xfem

# define constant averaging = "halfhalf"
define constant averaging = "hansbo"
#define constant averaging = "heaviside"

define constant lambda = 2
    
define constant heapsize = 1e8

# general constants
define constant zero = 0.0
define constant one = 1.0
define constant two = 2.0

define constant C1 = 1.0
define constant C2 = 4.0

#geometry constants
define constant r0 = 0.5+1e-12
define constant omega = 5
define constant omega2 = (1.2*pi) 

define coefficient lset
(
 sqrt(x*x+y*y)-(r0+0.2*sin(omega*atan2(x,y)))
)

numproc draw npd -coefficient=lset -label=levelset 

# henry weights
define constant bneg = 1.0
define constant bpos = 1.5

define constant aneg = 2.0
define constant apos = 1.0

define constant abneg = (aneg*bneg)
define constant abpos = (apos*bpos)

define coefficient solpos
(C1*sin(omega*atan2(x,y)+omega2)),
#((1-1.0/(sqrt(x*x+y*y)))*1.0*sin(omega*(atan2(x,y))+omega2*pi/180)),

define coefficient solneg
(C2*((x*x+y*y)-1.0/8.0)),

define coefficient rhspos
(bpos*C1*apos*omega*omega/(x*x+y*y)*sin(omega*atan2(x,y)+omega2)),

define coefficient rhsneg
(-bneg*4*C2*aneg),

define coefficient jumprhs 
(bpos*C1*sin(omega*atan2(x,y)+omega2)-bneg*C2*((x*x+y*y)-1.0/8.0)),

define coefficient fluxjumprhs
(
  (sqrt(x*x+y*y)/
   (
    sqrt(
         1+(omega/5.0*cos(omega*atan2(x,y)))/(r0+0.2*sin(omega*atan2(x,y)))*(omega/5.0*cos(omega*atan2(x,y)))/(r0+0.2*sin(omega*atan2(x,y)))
        )
   )
  )
  *
  (
   aneg*2*C2+apos*omega*C1*cos(omega*atan2(x,y)+omega2)*(omega/5.0*cos(omega*atan2(x,y)))/(r0+0.2*sin(omega*atan2(x,y)))*1.0/(x*x+y*y)
  )
)
,


define fespace fescomp
       -type=xstdfespace
       -type_std=h1ho 
       -order=1
       -dirichlet=[1,2]
       -ref_space=1
#       -dgjumps

numproc informxfem npix 
        -xstdfespace=fescomp
        -coef_levelset=lset

define gridfunction u -fespace=fescomp

define linearform f -fespace=fescomp # -print
xsource rhsneg rhspos
xnitscherhsjump_$(averaging) aneg apos bneg bpos jumprhs lambda
xnitscherhsfluxjump_$(averaging) bneg bpos fluxjumprhs
#xnitscherhsfluxjump_hansbo 1.0 1.0 1.0
# xsource solneg solpos

define bilinearform a -fespace=fescomp -printelmat #-eliminate_internal -keep_internal -symmetric -linearform=f # -printelmat -print
# xmass one one

xlaplace abneg abpos
xnitsche_$(averaging) aneg apos bneg bpos lambda

#xnitsche_minstab_hansbo aneg apos bneg bpos
#lo_ghostpenalty aneg apos one

numproc setvaluesx npsvx -gridfunction=u -coefficient_neg=solneg -coefficient_pos=solpos -boundary #-print

define preconditioner c -type=local -bilinearform=a #-test -block
#define preconditioner c -type=direct -bilinearform=a #-test
#define preconditioner c -type=bddc -bilinearform=a -test -block

##define preconditioner c -type=multigrid -bilinearform=a -test #-smoother=block

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6 # -print

define coefficient veczero
(0,0),

#numproc calccond npcc -bilinearform=a -inverse=pardiso #-symmetric

numproc xdifference npxd 
        -solution=u 
        -solution_n=solneg
        -solution_p=solpos
        -jumprhs=jumprhs
        -levelset=lset
        -interorder=5
        -henryweight_n=1
        -henryweight_p=1.5
        -diffusion_n=1
        -diffusion_p=1

numproc visualization npviz -scalarfunction=u -minval=-1.5 -maxval=1.5 #-comp=0

#numproc markinterface npmi -fespace=fescomp


coefficient err
(
  (sqrt(x*x+y*y)-(r0+0.2*sin(omega*atan2(x,y))) > 0) * (abs((u)-solpos)) 
 +(sqrt(x*x+y*y)-(r0+0.2*sin(omega*atan2(x,y))) < 0) * (abs((u)-solneg)) 
)

numproc draw npdraw -coefficient=err -label=error


