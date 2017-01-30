
# load geometry
geometry = d8_stokes.geo                                        
# and mesh
mesh = d8_stokes.vol.gz

#load xfem-library and python-bindings
shared = libngsxfem_xfem                                       
shared = libngsxfem_xstokes                                       
# pymodule = d7_stokes

define constant heapsize = 1e9

define constant R = 0.4
define constant one = 1.0

# interface description as zero-level
# define coefficient lset
# ( sqrt(x*x+y*y) - R),

define coefficient lset
( sqrt(x*x+y*y+2*z*z) - R),

# define coefficient lset
# ( x-0.127378 ),

define fespace fescomp_u                                            
       -type=xstdfespace                                      
       -type_std=h1ho                                          
       -order=2                 
       -dirichlet=[1]
       -empty
       # -dgjumps
       # -ref_space=5

numproc informxfem npi_uvx 
        -xstdfespace=fescomp_u
        -coef_levelset=lset


define fespace fescomp_p                                          
       -type=xstdfespace                                      
       -type_std=h1ho                                          
       -order=1                 
       # -empty
       # -ref_space=1                                           

numproc informxfem npi_px 
        -xstdfespace=fescomp_p
        -coef_levelset=lset

define fespace fescomp
       -type=compound
       -spaces=[fescomp_u,fescomp_u,fescomp_u,fescomp_p]

define gridfunction uvp -fespace=fescomp

define constant zero = 0.0
define constant one = 1.0
define constant lambda = 1000.0
define constant delta = 1.0

define coefficient s
0,1,0,0,

#numproc setvaluesx npsvx -gridfunction=uvp.2 -coefficient_neg=s -coefficient_pos=s -boundary

# integration on sub domains
define linearform f -fespace=fescomp
xsource one zero -comp=3
#xLBmeancurv one # naiv Laplace-Beltrami discretization 
#xmodLBmeancurv one lset # improved Laplace-Beltrami discretization 
# integration on sub domains
define bilinearform a -fespace=fescomp -symmetric -linearform=f -printelmat
xstokes one one 
# xlaplace one one -comp=1
# xnitsche one one one one lambda -comp=1
# xnitsche one one one one lambda -comp=2
# lo_ghostpenalty one one delta -comp=1
# lo_ghostpenalty one one delta -comp=2
# xmass one one -comp=1
# xmass 1.0 1.0 -comp=2

#define preconditioner c -type=local -bilinearform=a -test #-block           
define preconditioner c -type=direct -bilinearform=a -inverse=umfpack #-test 

numproc bvp npbvp -gridfunction=uvp -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6


define coefficient velocity ( (uvp.1.1, uvp.2.1, uvp.3.1) )
numproc draw npd1 -coefficient=velocity -label=velocity

define coefficient pressure ( uvp.4 )
numproc draw npd2 -coefficient=pressure -label=pressure

numproc draw npd3 -coefficient=lset -label=levelset

numproc visualization npviz 
        -scalarfunction=levelset
        -clipvec=[0,1,0]
        -centerpoint=[0,0,0]
#         # -vectorfunction=velocity
#         -minval=0 -maxval=0 
#         -nolineartexture -deformationscale=1 -subdivision=3


