
# load geometry
geometry = d7_stokes.in2d                                        
# and mesh
mesh = d7_stokes.vol.gz

#load xfem-library and python-bindings
shared = libngsxfem_xfem                                       
shared = libngsxfem_xstokes                                       
# pymodule = d7_stokes

define constant heapsize = 1e9

define constant R = 0.6666666
#define constant R = 0.5
define constant one = 1.0

define constant eps1 = 1e-1

# interface description as zero-level
define coefficient lset
( sqrt((x)*(x)+y*y) - R),

define fespace fescomp
       -type=xstokes
       -order=1               
       -dirichlet_vel=[1,2,3,4]
       #-empty_vel
       -dgjumps
       -ref_space=1

numproc informxstokes npi_px 
        -xstokesfespace=fescomp
        -coef_levelset=lset

define gridfunction uvp -fespace=fescomp
define gridfunction exu -fespace=fescomp
define gridfunction exuneg -fespace=fescomp

define constant mu1 = 1
define constant mu2 = 2

define coefficient alphapos
(1/mu2 + (1/mu1 - 1/mu2) * exp((x*x+y*y) - R*R)),

define coefficient alphaneg
(1/mu1),

define coefficient exactuxpos
(alphapos*exp(-1.0 * (x * x + y * y)) * -1.0 * y),

define coefficient exactuypos
(alphapos*exp(-1.0 *( x * x + y * y)) * x),

define coefficient exactuxneg
(alphaneg*exp(-1.0 * (x * x + y * y)) * -1.0 * y),

define coefficient exactuyneg
(alphaneg*exp(-1.0 *( x * x + y * y)) * x),

define coefficient exactp
(x*x*x),

define coefficient exactux
( (lset > 0) * exactuxpos + (lset < 0) * exactuxneg),

define coefficient exactuy
( (lset > 0) * exactuypos + (lset < 0) * exactuyneg),


numproc setvalues npsvex1 -gridfunction=exu.1.1 -coefficient=exactux
numproc setvalues npsvex2 -gridfunction=exu.2.1 -coefficient=exactuy
#numproc setvalues npsvex6 -gridfunction=exuneg.1.1 -coefficient=exactuxneg
#numproc setvalues npsvex7 -gridfunction=exuneg.2.1 -coefficient=exactuyneg
numproc setvalues npsvex5 -gridfunction=exu.3.1 -coefficient=exactp

numproc setvalues npsvex3 -gridfunction=uvp.1.1 -coefficient=exactux -boundary
numproc setvalues npsvex4 -gridfunction=uvp.2.1 -coefficient=exactuy -boundary

define constant zero = 0.0
define constant one = 1.0
#define constant none = -1.0
define constant lambda = 20
define constant delta = 1.0

define coefficient gammaf
0.0,
#(1.0/R),

define coefficient exactpneg
(x*x*x + (gammaf) - (pi*R*R/4.0*gammaf)),

define coefficient exactppos
(x*x*x - (pi*R*R/4.0*gammaf)),

define coefficient foneneg
(exp(-1* (x * x + y * y)) *((-8 * y) + (4 * x * x * y) + (4 * y * y * y))+ 3 * x * x),

define coefficient ftwoneg
(exp(-1* (x * x + y * y)) *((-4 * x * x * x) + (8 * x) - (4 * x * y * y))),

define coefficient fonepos
(exp(-1* (x * x + y * y)) * ((-8 * y) + (4 * x * x * y) + (4 * y * y * y))+ 3 * x * x),

define coefficient ftwopos
(exp(-1* (x * x + y * y)) * ((-4 * x * x * x) + (8 * x) - (4 * x * y * y))),

define coefficient ghost
-0.1,

define coefficient eps
(1e-6),

#numproc setvaluesx npsvx -gridfunction=uvp.2 -coefficient_neg=s -coefficient_pos=s -boundary

# integration on sub domains
define linearform f -fespace=fescomp
xsource foneneg fonepos -comp=1
xsource ftwoneg ftwopos -comp=2

#xsource zero zero -comp=2
#xGammaForce gammaf
#xLBmeancurv one # naiv Laplace-Beltrami discretization 
#xmodLBmeancurv one lset # improved Laplace-Beltrami discretization 
# integration on sub domains
define bilinearform a -fespace=fescomp -linearform=f -printelmat
xstokes mu1 mu2
myghostpenalty ghost -comp=3
xstokesnitsche mu1 mu2 lambda

#xmass one one -comp=1
#xmass one one -comp=2

#myghostpenalty ghost -comp=2
#myghostpenalty ghost -comp=1
#xmass eps eps -comp=3
# xlaplace one one -comp=1
# xnitsche one one one one lambda -comp=1
# xnitsche one one one one lambda -comp=2
#lo_ghostpenalty one one ghost -comp=3
# lo_ghostpenalty one one delta -comp=2


#define preconditioner c -type=local -bilinearform=a -test #-block           
define preconditioner c -type=direct -bilinearform=a -inverse=pardiso #-test 

numproc bvp npbvp -gridfunction=uvp -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6


define coefficient velocity ( (uvp.1, uvp.2) )
numproc draw npd1 -coefficient=velocity -label=velocity

define coefficient pressure ( (uvp.3) )
numproc draw npd2 -coefficient=pressure -label=pressure

numproc draw npd3 -coefficient=lset -label=levelset

numproc visualization npviz 
        -scalarfunction=levelset
        # -vectorfunction=velocity
        -minval=0 -maxval=0 
        -nolineartexture -deformationscale=1 -subdivision=3

define bilinearform b1 -fespace=fescomp -symmetric -nonassemble
xvis one -comp=2
#xmass one one -comp=2
#xmass one one -comp=3

define fespace feerror -type=l2ho -order=0

define gridfunction erroru -fespace=feerror

#numproc drawflux test -bilinearform=b1 -solution=uvp -label=qqqq

numproc difference checkdiff -bilinearform1=b1 -solution1=uvp -bilinearform2=b1 -solution2=exu -diff=erroru

#numproc difference checkdiff -bilinearform=b1 -solution=uvp -function=fone -diff=erroru



# numproc xdifference diffp -solution=uvp.3 -solution_n=exactpneg -solution_p=exactppos
numproc xdifference diffp2 -solution=uvp.3 -solution_n=exactpneg -solution_p=exactppos -nooutput
numproc xdifference diffp3 -solution=uvp.3 -solution_n=exactpneg -solution_p=exactppos

numproc xdifference diffpuv1 -solution=uvp.1 -solution_n=exactuxneg -solution_p=exactuxpos
numproc xdifference diffpuv2 -solution=uvp.2 -solution_n=exactuyneg -solution_p=exactuypos
        
numproc evaluate eval1 -bilinearform=b1 -gridfunction=uvp -point=[-1,0,0] -point2=[0,0,0] -filename=evalplot
numproc evaluate eval2 -bilinearform=b1 -gridfunction=exu -point=[-1,0,0] -point2=[0,0,0] -filename=evalplot2
#-reference=exu.3


define fespace fesvelstd
        -type=h1ho
        -order=2
define fespace fespstd
        -type=h1ho
        -order=1

define fespace fesvelnegpos
        -type=compound
        -spaces=[fesvelstd,fesvelstd]

define fespace fespnegpos
        -type=compound
        -spaces=[fespstd,fespstd]

define gridfunction gf_u_negpos -fespace=fesvelnegpos -novisual
define gridfunction gf_v_negpos -fespace=fesvelnegpos -novisual
define gridfunction gf_p_negpos -fespace=fespnegpos -novisual

numproc xtonegpos npxtonegposu -xstd_gridfunction=uvp.1 -negpos_gridfunction=gf_u_negpos
numproc xtonegpos npxtonegposv -xstd_gridfunction=uvp.2 -negpos_gridfunction=gf_v_negpos
numproc xtonegpos npxtonegposp -xstd_gridfunction=uvp.3 -negpos_gridfunction=gf_p_negpos

numproc vtkoutput npout -filename=task4
        -coefficients=[lset]
        -gridfunctions=[gf_u_negpos.1,gf_u_negpos.2,gf_v_negpos.1,gf_v_negpos.2,gf_p_negpos.1,gf_p_negpos.2]
        -fieldnames=[levelset,uneg,upos,vneg,vpos,pneg,ppos]        