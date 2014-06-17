
# load geometry
geometry = d1_simple.in2d

# and mesh
mesh = d4_reg8.vol.gz
#mesh = d1b_reg.vol.gz
# mesh = d1_simple.vol.gz


shared = libngsxfem_xfem

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
define constant omega2 = (1*pi) 

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

define constant lambda = 2

define fespace fescomp
       -type=xh1fespace
       -order=1
       -dirichlet=[1,2]
       -ref_space=1
#       -dgjumps

numproc informxfem npix 
        -xh1fespace=fescomp
        -coef_levelset=lset

define gridfunction u -fespace=fescomp

define linearform f -fespace=fescomp # -print
xsource rhsneg rhspos
xnitscherhsjump_hansbo aneg apos bneg bpos jumprhs lambda
xnitscherhsfluxjump_hansbo bneg bpos fluxjumprhs

# xsource solneg solpos

define bilinearform a -fespace=fescomp -printelmat #-eliminate_internal -keep_internal -symmetric -linearform=f # -printelmat -print
# xmass one one

xlaplace abneg abpos
xnitsche_hansbo aneg apos bneg bpos lambda

#xnitsche_minstab_hansbo aneg apos bneg bpos
#lo_ghostpenalty aneg apos one

numproc setvaluesx npsvx -gridfunction=u -coefficient_neg=solneg -coefficient_pos=solpos -boundary -print

#define preconditioner c -type=local -bilinearform=a -test -block
define preconditioner c -type=direct -bilinearform=a #-test
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

numproc visualization npviz -scalarfunction=u #-comp=0

# numproc xoutput npxo 
#         -solution=u 
#         -solution_n=solneg 
#         -solution_p=solpos 
#         -subdivision=1
# #        -showerror
#         -overlapeps=0e-2

#         -negcolor=green!40 #2 #white!20!black
#         -fillnegopacity=1.0

#         -edgenegcolor=black #2
#         -edgenegstyle=thick
#         -drawnegopacity=1.0

#         -fineedgenegcolor=white!20!black #2
#         -fineedgenegstyle=dashed
#         -fineedgenegopacity=0.3

#         -poscolor=white!40!black #1
#         -fillposopacity=0.1

#         -edgeposcolor=black #1
#         -edgeposstyle=thick
#         -drawposopacity=1.0

#         -fineedgeposcolor=black #1
#         -fineedgeposstyle=dashed
#         -fineedgeposopacity=0.3

#         -viewpoint=(6.0,9.0,-3.0)
#         -lookat=(1.0,0.0,1.0)
#         -scalex=10.0        
#         -scaley=10.0        
#         -scalez=2.0

#         # -viewpoint=(0.0,5.0,5.0)
#         # -lookat=(0.0,0.0,0.0)
#         # -scalex=20.0        
#         # -scaley=20.0        
#         # -scalez=0.0
