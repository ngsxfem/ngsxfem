
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

#geometry constants
define constant r0 = 0.5
define constant omega = 5

define coefficient lset
(
  sqrt(x*x+y*y)-(r0+0.2*sin(omega*atan2(x,y)))
)

numproc draw npd -coefficient=lset -label=levelset 

# henry weights
define constant bneg = 1.0
define constant bpos = 1.5

define constant aneg = 1.0
define constant apos = 1.5

define constant abneg = (aneg*bneg)
define constant abpos = (apos*bpos)

define coefficient f
(1-0.25*(x*x+y*y)),

# define coefficient df
# (cos(x+y)),

define coefficient mddf
1.0,

define coefficient rhspos
0, #(bpos*apos),

define coefficient rhsneg
1, #(bpos*aneg),

define coefficient solpos
0, #(1-0.25*(x*x+y*y)),

# define coefficient solposx
# ((df)),

# define coefficient solposy
# ((df)),

define coefficient solneg
(bpos/bneg*(1-0.25*(x*x+y*y))),

# define coefficient solnegx
# (bpos/bneg*(df)),

# define coefficient solnegy
# (bpos/bneg*(df)),

define constant lambda = 2

define fespace fescomp
       -type=xh1fespace
       -order=2
       -dirichlet=[1,2]
       -ref_space=2
#       -dgjumps

numproc informxfem npix 
        -xh1fespace=fescomp
        -coef_levelset=lset

define gridfunction u -fespace=fescomp

define linearform f -fespace=fescomp # -print
xsource rhsneg rhspos
#xsource solneg solpos

define bilinearform a -fespace=fescomp -printelmat #-eliminate_internal -keep_internal -symmetric -linearform=f # -printelmat -print
#xmass one one
xlaplace abneg abpos
#xnitsche_minstab_hansbo aneg apos bneg bpos
xnitsche_hansbo aneg apos bneg bpos lambda
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


############### REFERENCE SOLUTION ###########
define fespace fescomp2
       -type=xh1fespace
       -order=4
       -dirichlet=[1,2]
       -ref_space=2
#       -dgjumps

numproc informxfem npix2
        -xh1fespace=fescomp2
        -coef_levelset=lset

define gridfunction u2 -fespace=fescomp2

define linearform f2 -fespace=fescomp2 # -print
xsource rhsneg rhspos
#xsource solneg solpos

define bilinearform a2 -fespace=fescomp2 #-printelmat #-eliminate_internal -keep_internal -symmetric -linearform=f # -printelmat -print
#xmass one one
xlaplace abneg abpos
#xnitsche_minstab_hansbo aneg apos bneg bpos
xnitsche_hansbo aneg apos bneg bpos lambda
#lo_ghostpenalty aneg apos one

numproc setvaluesx npsvx2 -gridfunction=u2 -coefficient_neg=solneg -coefficient_pos=solpos -boundary #-print

#define preconditioner c2 -type=local -bilinearform=a -test #-block
define preconditioner c2 -type=direct -bilinearform=a2 #-test
#define preconditioner c2 -type=bddc -bilinearform=a #-test -block

##define preconditioner c -type=multigrid -bilinearform=a -test #-smoother=block

numproc bvp npbvp2 -gridfunction=u2 -bilinearform=a2 -linearform=f2 -solver=cg -preconditioner=c2 -maxsteps=1000 -prec=1e-6 # -print


################################


numproc xdifference npxd 
        -solution=u 
        -solution2=u2
        -solution_n=solneg
        -solution_p=solpos
        -levelset=lset
        -interorder=5
        -henryweight_n=1
        -henryweight_p=1.5

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
