# Example:
# Streamline Diffusion combined with Nitsche-XFEM
# ... artificial interface
#
# The interface is prescribed with the coef "lset",
# ...
#
# Solves the problem:
#
# ...TODO...
#
# Details:
#
# → ....
# → ....
#
# Things to try here:
#   1. ...
#   2. ... switch scalings ...
#

geometry = d3_sdxnitsche.in2d
mesh = d3_sdxnitsche.vol

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
       -type=xstdfespace
       -type_std=h1ho
       -order=1
#       -dirichlet=[4]
       -ref_space=0
#       -dgjumps

numproc informxfem npix 
        -xstdfespace=vh1x
        -coef_levelset=lset


define coefficient one
1,

define gridfunction uh1x -fespace=vh1x

define constant a_n = (9e-7)
define constant a_p = (4e-7)
define constant b_n = 27.0
define constant b_p = 11.0

define constant gamma_tr = 0.01

define constant dn = 0.2
define constant Bn = (-16.0/27.0)    # delta A
define constant Cn = (-log(gamma_tr)/dn)

define constant dp = 0.15
define constant Bp = 1.0
define constant Cp = (-log(gamma_tr)/dp)

define constant bndpen = 1e9

define coefficient brhsneg
( b_n * ( (Bn*Cn/x* (-0.5*y/sqrt(x) - a_n*(Cn-0.25*y/x*(-Cn*y/x-3/sqrt(x))) )) * exp(Cn/sqrt(x)*y)) ),

define coefficient brhspos
( b_p * ( (Bp*Cp/x* (0.5*y/sqrt(x) - a_p*(Cp+0.25*y/x*(Cp*y/x-3/sqrt(x))) )) * exp(-Cp/sqrt(x)*y)) ),

define coefficient rhsneg
( ( (Bn*Cn/x* (-0.5*y/sqrt(x) - a_n*(Cn-0.25*y/x*(-Cn*y/x-3/sqrt(x))) )) * exp(Cn/sqrt(x)*y)) ),

define coefficient rhspos
( ( (Bp*Cp/x* (0.5*y/sqrt(x) - a_p*(Cp+0.25*y/x*(Cp*y/x-3/sqrt(x))) )) * exp(-Cp/sqrt(x)*y)) ),

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

define coefficient w_pos
(1,0),

define coefficient w_neg
(1,0),

define coefficient bw_pos
b_p,0,

define coefficient bw_neg
b_n,0,

define coefficient rob_n
0,0,0,(bndpen),

define coefficient rob_p
0,0,0,(bndpen),

define coefficient sol_n
(1.0-16.0/27.0*exp(Cn/sqrt(x)*y)),

define coefficient solx_n
(-16.0/27.0*Cn/sqrt(x)*exp(Cn/sqrt(x)*y)*(-0.5*y/x)),

define coefficient soly_n
(-16.0/27.0*Cn/sqrt(x)*exp(Cn/sqrt(x)*y)*(1.0)),

define coefficient sol_p
(exp(-Cp/sqrt(x)*y)),

define coefficient solx_p
(-1.0*Cp/sqrt(x)*exp(-Cp/sqrt(x)*y)*(-0.5*y/x)),

define coefficient soly_p
(-1.0*Cp/sqrt(x)*exp(-Cp/sqrt(x)*y)*(1.0)),

define coefficient neu_n
(-a_n*16.0/27.0*Cn/sqrt(x)*exp(Cn/sqrt(x)*y)),
(a_n*(8.0/27.0)*Cn/sqrt(x)*exp(Cn/sqrt(x)*y)*(y/x)),
(a_n*16.0/27.0*Cn/sqrt(x)*exp(Cn/sqrt(x)*y)),
(bndpen * (1.0-16.0/27.0*exp(Cn/sqrt(x)*y))),

define coefficient neu_p
(-a_p*Cp/sqrt(x)*exp(-Cp/sqrt(x)*y)),
( a_p*Cp/sqrt(x)*exp(-Cp/sqrt(x)*y)*(0.5*y/x)),
( a_p*Cp/sqrt(x)*exp(-Cp/sqrt(x)*y)),
(bndpen * exp(-Cp/sqrt(x)*y)),

numproc setvaluesx npsvx -gridfunction=uh1x -coefficient_neg=sol_n -coefficient_pos=sol_p -boundary #-print

define bilinearform m -fespace=vh1x # -symmetric
xmass b_n b_p 

define coefficient small
1e-2,

define bilinearform a -fespace=vh1x #-symmetric
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
xrobin rob_n rob_p


#lo_ghostpenalty small

define linearform f -fespace=vh1x
xneumann neu_n neu_p
xsource brhsneg brhspos
sdxsource betaneg betapos 
          alphaneg alphapos 
          w_neg w_pos 
          zero zero
          rhsneg rhspos 

#define preconditioner cb -type=bddc -bilinearform=a -test
#define preconditioner c -type=local -bilinearform=a #-test
define preconditioner c -type=direct -bilinearform=a -inverse=pardiso

numproc bvp npx -bilinearform=a -linearform=f -gridfunction=uh1x -solver=cg -preconditioner=c -prec=1e-10 -maxsteps=5000 -print

numproc xdifference npxd 
        -solution=uh1x 
        -solution_n=sol_n 
        -solution_p=sol_p 
        -levelset=lset 
        -intorder=4 
        -threshold=0.1
        -henryweight_p=11.0 
        -henryweight_n=27.0
        -diffusion_n=9e-7
        -diffusion_p=4e-7
        -convection_n=w_neg
        -convection_p=w_pos


numproc visualization npviz -scalarfunction=uh1x #-comp=0


numproc xoutput npxo 
        -solution=uh1x
        -solution_n=sol_n
        -solution_p=sol_p 
        -subdivision=1
#        -showerror
        -overlapeps=0e-2

        -negcolor=green!40 #2 #white!20!black
        -fillnegopacity=1.0

        -edgenegcolor=black #2
        -edgenegstyle=thick
        -drawnegopacity=1.0

        -fineedgenegcolor=white!20!black #2
        -fineedgenegstyle=dashed
        -fineedgenegopacity=0.3

        -poscolor=white!40!black #1
        -fillposopacity=0.1

        -edgeposcolor=black #1
        -edgeposstyle=thick
        -drawposopacity=1.0

        -fineedgeposcolor=black #1
        -fineedgeposstyle=dashed
        -fineedgeposopacity=0.3

        -viewpoint=(6.0,9.0,-3.0)
        -lookat=(1.0,0.0,1.0)
        -scalex=10.0        
        -scaley=10.0        
        -scalez=2.0

        # -viewpoint=(0.0,5.0,5.0)
        # -lookat=(0.0,0.0,0.0)
        # -scalex=20.0        
        # -scaley=20.0        
        # -scalez=0.0
