
# load geometry
geometry = d7_stokes.in2d                                        
# and mesh
mesh = d7_stokes.vol.gz

#load xfem-library and python-bindings
shared = libngsxfem_xfem                                       
shared = libngsxfem_xstokes                                       
# pymodule = d7_stokes
flags tracer = -max_size=0
define constant heapsize = 1e9

define constant R = 0.666666666666666
#define constant R = 0.5
define constant one = 1.0

define constant eps1 = 1e-1

# interface description as zero-level
define coefficient lset
( sqrt((x)*(x)+y*y) - R),

define fespace fescomp
       -type=xstokes
       -order=2
       -dirichlet_vel=[1,2,3,4]
       #-empty_vel
       -dgjumps
       -ref_space=3

numproc informxstokes npi_px 
        -xstokesfespace=fescomp
        -coef_levelset=lset

define gridfunction uvp -fespace=fescomp
define gridfunction exu -fespace=fescomp
define gridfunction exuneg -fespace=fescomp

define constant mu1 = 1
define constant mu2 = 100


define constant zero = 0.0
define constant one = 1.0
#define constant none = -1.0
define constant lambda = 20
define constant delta = 1.0

define coefficient gammaf
0.0,
#(1.0/R),

define coefficient foneneg
-1.0,

define coefficient ftwoneg
0.0,

define coefficient ghost
-0.01,

define coefficient eps
(1e-6),

#numproc setvaluesx npsvx -gridfunction=uvp.2 -coefficient_neg=s -coefficient_pos=s -boundary

# integration on sub domains
define linearform f -fespace=fescomp
xsource ftwoneg foneneg -comp=2


define bilinearform a -fespace=fescomp -linearform=f -printelmat
xstokes mu1 mu2
xstokesnitsche mu1 mu2 lambda
lo_ghostpenalty one one ghost -comp=3


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

#numproc xdifference diffp2 -solution=uvp.3 -solution_n=zero -solution_p=zero -nooutput

        
numproc xtonegpos npxtonegposu -xstd_gridfunction=uvp.1 -negpos_gridfunction=gf_u_negpos
numproc xtonegpos npxtonegposv -xstd_gridfunction=uvp.2 -negpos_gridfunction=gf_v_negpos
numproc xtonegpos npxtonegposp -xstd_gridfunction=uvp.3 -negpos_gridfunction=gf_p_negpos

numproc vtkoutput npout -filename=bubble2_
        -coefficients=[lset]
        -gridfunctions=[gf_u_negpos.1,gf_u_negpos.2,gf_v_negpos.1,gf_v_negpos.2,gf_p_negpos.1,gf_p_negpos.2]
        -fieldnames=[levelset,uneg,upos,vneg,vpos,pneg,ppos]
        -subdivision=1

          