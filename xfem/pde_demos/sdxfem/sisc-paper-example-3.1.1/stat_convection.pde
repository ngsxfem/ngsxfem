geometry = square_conv.in2d
mesh = square_conv_240els.vol.gz

shared = libngsxfem_xfem

define constant heapsize = 1e9

# general constants
define constant zero = 0.0
define constant one = 1.0
define constant two = 2.0

#geometry constants
define constant r0 = 0.5
define constant omega = 5

define coefficient lset
(y),

numproc draw npd -coefficient=lset -label=levelset 

define fespace vh1x
       -type=xh1fespace
       -order=1
       -dirichlet=[1,2,3,4]
       -ref_space=0
#       -dgjumps

numproc informxfem npix 
        -xh1fespace=vh1x
        -coef_levelset=lset


define coefficient one
1,

define gridfunction uh1x -fespace=vh1x

define constant a_n = (2e-7)
define constant a_p = (1e-7)
define constant b_n = 3.0
define constant b_p = 2.0

define coefficient one
1,

define coefficient rhsneg
( ( pi*pi*4/3*sin(pi*(x+y))*a_n + 2*pi/3*cos(pi*(x+y)) ) ),

define coefficient rhspos
( ( pi*pi*(25/9)*sin(pi*(x+4/3*y))*a_p + pi * cos(pi*(x+4/3*y)) ) ),

define coefficient brhsneg
( b_n * ( pi*pi*4/3*sin(pi*(x+y))*a_n + 2*pi/3*cos(pi*(x+y)) ) ),

define coefficient brhspos
( b_p * ( pi*pi*(25/9)*sin(pi*(x+4/3*y))*a_p + pi * cos(pi*(x+4/3*y)) ) ),

define coefficient alphaneg
(a_n),

define coefficient alphapos
(a_p),

define coefficient alphabetaneg
(a_n*b_n),

define coefficient alphabetapos
(a_p*b_p),

define coefficient betaneg
(b_n),

define coefficient betapos
(b_p),

define coefficient betamassneg
(b_n*1),

define coefficient betamasspos
(b_p*1),

define coefficient lambda
5,

define coefficient zero
0,

define coefficient bw_pos
(b_p,0),

define coefficient bw_neg
(b_n,0),

define coefficient w_pos
(1,0),

define coefficient w_neg
(1,0),

define coefficient sol_n
(2/3*sin(pi*(x+y)) ),

define coefficient sol_p
(sin(pi*(x+4/3*y)) ),

define constant bndpen = 1e9

define coefficient neu_sol_n
(bndpen*(sol_n) ),
(bndpen*(sol_n) ),
(bndpen*(sol_n) ),
(bndpen*(sol_n) ),

define coefficient neu_sol_p
(bndpen*(sol_p) ),
(bndpen*(sol_p) ),
(bndpen*(sol_p) ),
(bndpen*(sol_p) ),

define coefficient rob_n
(bndpen),(bndpen),(bndpen),(bndpen),

define coefficient rob_p
(bndpen),(bndpen),(bndpen),(bndpen),

numproc setvaluesx npsvx -gridfunction=uh1x -coefficient_neg=sol_n -coefficient_pos=sol_p -boundary #-print

define bilinearform m -fespace=vh1x # -symmetric
xmass b_n b_p 


define coefficient small
1e-2,

define bilinearform a -fespace=vh1x # -symmetric
xrobin rob_n rob_p
xlaplace alphabetaneg alphabetapos 
xnitsche_conv alphaneg alphapos 
              betaneg betapos 
              lambda 
              w_neg w_pos


xconvection bw_neg bw_pos

sdstab betaneg betapos 
       alphaneg alphapos 
       w_neg w_pos 
       zero zero

#lo_ghostpenalty small


define linearform f -fespace=vh1x
xneumann neu_sol_n neu_sol_p

xsource brhsneg brhspos

sdxsource betaneg betapos 
          alphaneg alphapos 
          w_neg w_pos 
          zero zero
          rhsneg rhspos 

define linearform f2 -fespace=vh1x
xsource rhsneg rhspos

#define preconditioner c -type=bddc -bilinearform=a #-block #-test
define preconditioner c -type=local -bilinearform=a #-block #-test
#define preconditioner c -type=direct -bilinearform=a -inverse=pardiso
#define preconditioner c -type=direct -bilinearform=a -inverse=sparsecholesky

numproc bvp npx -bilinearform=a -linearform=f -gridfunction=uh1x -solver=gmres -preconditioner=c -prec=1e-10 -maxsteps=5000 -print

#numproc bvp npx -bilinearform=m -linearform=f2 -gridfunction=uh1x -preconditioner=c

numproc xdifference npxd 
        -solution=uh1x
        -solution_n=sol_n
        -solution_p=sol_p
        -levelset=lset
        -intorder=4
        #-threshold=-0.1
        -henryweight_n=3
        -henryweight_p=2
        -diffusion_n=2e-7
        -diffusion_p=1e-7
        -convection_n=w_neg
        -convection_p=w_pos


numproc visualization npviz -scalarfunction=uh1x #-comp=0