
# load geometry
geometry = d7_stokes.in2d                                        
# and mesh
mesh = d7_stokes.vol.gz

#load xfem-library and python-bindings
shared = libngsxfem_xfem                                       
shared = libngsxfem_xstokes                                       
# pymodule = d7_stokes

define constant heapsize = 1e9

define constant R = 0.75001
define constant one = 1.0

define constant eps1 = 1e-1
define constant R = 0.5+eps1

# interface description as zero-level
define coefficient lset
( sqrt((x)*(x)+y*y) - R),

define fespace fescomp
       -type=xstokes
       -order=1                 
       -dirichlet_vel=[1,2,3,4]
       -empty_vel
       -dgjumps
       -ref_space=1

numproc informxstokes npi_px 
        -xstokesfespace=fescomp
        -coef_levelset=lset

define gridfunction uvp -fespace=fescomp
define gridfunction exu -fespace=fescomp

define coefficient exactuxpos
(exp(-1.0 * (x * x + y * y)) * -1.0 * y ),

define coefficient exactuypos
(exp(-1.0 *( x * x + y * y)) * x ),

define coefficient exactp
(x*x*x),


numproc setvalues npsvex1 -gridfunction=exu.1.1 -coefficient=exactuxpos
numproc setvalues npsvex2 -gridfunction=exu.2.1 -coefficient=exactuypos
numproc setvalues npsvex5 -gridfunction=exu.3.1 -coefficient=exactp

numproc setvalues npsvex3 -gridfunction=uvp.1.1 -coefficient=exactuxpos -boundary
numproc setvalues npsvex4 -gridfunction=uvp.2.1 -coefficient=exactuypos -boundary

define constant zero = 0.0
define constant one = 1.0
#define constant none = -1.0
define constant lambda = 1000.0
define constant delta = 1.0
define constant mu1 = 1
define constant mu2 = 1

define coefficient gammaf
(0.0),
#(1.0/R),

define coefficient exactpneg
(x*x*x + (gammaf) - (pi*R*R/4.0*gammaf)),

define coefficient exactppos
(x*x*x - (pi*R*R/4.0*gammaf)),

define coefficient fone
(exp(-1* (x * x + y * y))  * ((-8 * y) + (4 * x * x * y) + (4 * y * y * y))+ 3 * x * x),

define coefficient ftwo
(exp(-1* (x * x + y * y)) * ((-4 * x * x * x) + (8 * x) - (4 * x * y * y))),

define coefficient ghost
-0.01,

define coefficient eps
(1e-6),

#numproc setvaluesx npsvx -gridfunction=uvp.2 -coefficient_neg=s -coefficient_pos=s -boundary

# integration on sub domains
define linearform f -fespace=fescomp
xsource fone fone -comp=1
xsource ftwo ftwo -comp=2

#xsource zero zero -comp=2
xGammaForce gammaf
#xLBmeancurv one # naiv Laplace-Beltrami discretization 
#xmodLBmeancurv one lset # improved Laplace-Beltrami discretization 
# integration on sub domains
define bilinearform a -fespace=fescomp -symmetric -linearform=f -printelmat
xstokes mu1 mu2 
myghostpenalty ghost -comp=3
#xmass eps eps -comp=3
# xlaplace one one -comp=1
# xnitsche one one one one lambda -comp=1
# xnitsche one one one one lambda -comp=2
#lo_ghostpenalty one one ghost -comp=3
# lo_ghostpenalty one one delta -comp=2
# xmass one one -comp=1
# xmass 1.0 1.0 -comp=2

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
xvis one -comp=1
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

#-reference=exu.3
