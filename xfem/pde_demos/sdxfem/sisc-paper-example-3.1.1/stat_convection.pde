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
#       -dirichlet=[1,2,3,4]
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

define constant mm = 0.0 #mass term?!
define constant mb_n = (b_n*mm)
define constant mb_p = (b_p*mm)

define constant pen = 1e9

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

define coefficient pen
0,0,0,0,

define coefficient cpen
0,0,0,0,

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

define coefficient solx_n
(2/3*pi*cos(pi*(x+y)) ),

define coefficient soly_n
(2/3*pi*cos(pi*(x+y)) ),

define coefficient sol_p
(sin(pi*(x+4/3*y)) ),

define coefficient solx_p
(    pi*cos(pi*(x+4/3*y)) ),

define coefficient soly_p
(4/3*pi*cos(pi*(x+4/3*y)) ),

define coefficient neu_sol_n
(2/3*pen*sin(pi*(x+y)) ),
(2/3*pen*sin(pi*(x+y)) ),
(2/3*pen*sin(pi*(x+y)) ),
(2/3*pen*sin(pi*(x+y)) ),

define coefficient neu_sol_p
(pen*sin(pi*(x+4/3*y)) ),
(pen*sin(pi*(x+4/3*y)) ),
(pen*sin(pi*(x+4/3*y)) ),
(pen*sin(pi*(x+4/3*y)) ),

define coefficient rob_n
pen,
pen,
pen,
pen,

define coefficient rob_p
pen,
pen,
pen,
pen,




numproc setvaluesx npsvx -gridfunction=uh1x -coefficient_neg=sol_n -coefficient_pos=sol_p -boundary #-print

define bilinearform m -fespace=vh1x # -symmetric
xmass b_n b_p 


define coefficient small
1e-2,

define bilinearform a -fespace=vh1x # -symmetric
xrobin rob_n rob_p
xlaplace alphabetaneg alphabetapos 
xnitsche_hansbo alphaneg alphapos 
                betaneg betapos 
                lambda

#xnitscheconv alphaneg alphapos 
#             betaneg betapos 
#             w_neg w_pos
#             lambda

xconvection bw_neg bw_pos

sdstab betaneg betapos 
       alphaneg alphapos 
       w_neg w_pos 
       mb_n mb_p

xmass mb_n mb_p

#lo_ghostpenalty small


define linearform f -fespace=vh1x
xneumann neu_sol_n neu_sol_p

xsource brhsneg brhspos

sdxsource betaneg betapos 
          alphaneg alphapos 
          w_neg w_pos 
          mb_n mb_p
          rhsneg rhspos 

define linearform f2 -fespace=vh1x
xsource rhsneg rhspos

#define preconditioner cb -type=bddc -bilinearform=a -test
#define preconditioner cl -type=local -bilinearform=a -test
define preconditioner c -type=direct -bilinearform=a -inverse=pardiso

numproc bvp npx -bilinearform=a -linearform=f -gridfunction=uh1x -preconditioner=c -prec=1e-22

#numproc bvp npx -bilinearform=m -linearform=f2 -gridfunction=uh1x -preconditioner=c

numproc xdifference npxd 
        -solution=uh1x
        -solution_n=sol_n
        -solution_p=sol_p
        -levelset=lset
        -interorder=4
        #-threshold=-0.1
        -henryweight_n=3
        -henryweight_p=2
        -diffusion_n=2e-7
        -diffusion_p=1e-7

#numproc xdifference npxd 
#        -solution=uh1x 
#        -function_n=sol_n 
#        -function_p=sol_p 
#        -derivative_x_n=solx_n -derivative_y_n=soly_n 
#        -derivative_x_p=solx_p -derivative_y_p=soly_p
#        -levelset=lset -intorder=4 -threshold=-0.1
#        -henryweight_p=2 - henryweight_n=3


#define bilinearform mh1 -fespace=vh1x -nonassemble
#mass one -comp=1

#define bilinearform mx -fespace=vh1x -nonassemble
#mass one -comp=2

#define bilinearform mvish1 -fespace=vh1x -nonassemble
#visx one 

#define bilinearform mvish1_lset -fespace=vh1x -nonassemble
#visx_lset lset 

#define bilinearform mvish1_lset_cut -fespace=vh1x -nonassemble
#visx_lset_sign lset 

#numproc drawflux npdf0 -bilinearform=mvish1 -solution=uh1x -label=u
#numproc drawflux npdf0 -bilinearform=mvish1_lset -solution=uh1x -label=u  -applyd
#numproc drawflux npdf1 -bilinearform=mh1 -solution=uh1x -label=u_h1
#numproc drawflux npdf2 -bilinearform=mx -solution=uh1x -label=u_x
#numproc drawflux npdf3 -bilinearform=mvish1_lset -solution=uh1x -label=u_part
#numproc drawflux npdf4 -bilinearform=mvish1_lset_cut -solution=uh1x -label=u_part_cut -applyd

#numproc visualization mpviz -scalarfunction=u -comp=1 -deformationscale=1 -subdivision=2

numproc visualization npviz -scalarfunction=uh1x #-comp=0