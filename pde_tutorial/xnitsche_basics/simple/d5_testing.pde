
# load geometry
geometry = d5_simple.in2d

# and mesh
mesh = d5_simple.vol.gz

shared = libngsxfem_xfem

define constant heapsize = 1e8

# general constants
define constant zero = 0.0
define constant one = 1.0
define constant two = 2.0

#geometry constants
define constant r0 = 0.33333
define constant omega = 5
define constant omega2 = 18  # winkel in grad

# define constant omega = 0
# define constant omega2 = 90 # winkel in grad

# define coefficient lset
# (
#   sqrt(x*x+y*y)-(r0+0.2*sin(omega*atan2(x,y)))
# )

define coefficient lset
(
  sqrt(x*x)-r0
)


numproc draw npd -coefficient=lset -label=levelset 

# henry weights
define constant bneg = 1.0
define constant bpos = 1.0

define constant aneg = 1.0
define constant apos = 1.0

define constant abneg = (aneg*bneg)
define constant abpos = (apos*bpos)

define coefficient solpos
(abs(x)),
#((1-1.0/(sqrt(x*x+y*y)))*1.0*sin(omega*(atan2(x,y))+omega2*pi/180)),

# define coefficient solposx
# ((df)),

# define coefficient solposy
# ((df)),

define coefficient solneg
r0,
#(x*x+y*y),

# define coefficient solnegx
# (bpos/bneg*(df)),

# define coefficient solnegy
# (bpos/bneg*(df)),

define coefficient rhspos
0,

#*(omega*omega*(1.0/(x*x+y*y)-1.0/((x*x+y*y)*sqrt(x*x+y*y))) - (1.0/(sqrt(x*x+y*y))-1.0/((x*x+y*y)*sqrt(x*x+y*y))))

define coefficient rhsneg
0,
#(-4*aneg),

define coefficient jumprhs 
0,
#(bpos*sin(omega*atan2(x,y)+omega2*pi/180)),#-bneg*(x*x+y*y)),

define coefficient fluxjumprhs
-1,
#(2*aneg*r0), 



define constant lambda = 2

define fespace fescomp
       -type=xh1fespace
       -order=1
       -dirichlet=[2]
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
        -henryweight_p=1
        -diffusion_n=1
        -diffusion_p=1

numproc visualization npviz -scalarfunction=u #-comp=0

# numproc xoutput npxo 
#         -solution=u 
#         -solution_n=solneg 
#         -solution_p=solpos 
#         -subdivision=2
# #        -showerror
#         -overlapeps=0e-2

#         -negcolor=white!80!black2 #white!20!black
#         -fillnegopacity=1.0

#         -edgenegcolor=black2
#         -edgenegstyle=thick
#         -drawnegopacity=1.0

#         -fineedgenegcolor=white!20!black2
#         -fineedgenegstyle=dashed
#         -fineedgenegopacity=0.5

#         -poscolor=white!40!black1
#         -fillposopacity=1.0

#         -edgeposcolor=black1
#         -edgeposstyle=thick
#         -drawposopacity=1.0

#         -fineedgeposcolor=black1
#         -fineedgeposstyle=dashed
#         -fineedgeposopacity=0.5

#         -viewpoint=(6.0,6.0,-3.0)
#         -lookat=(1.0,0.0,1.0)
#         -scalex=10.0        
#         -scaley=10.0        
#         -scalez=10.0

#         # -viewpoint=(0.0,5.0,5.0)
#         # -lookat=(0.0,0.0,0.0)
#         # -scalex=20.0        
#         # -scaley=20.0        
#         # -scalez=0.0
