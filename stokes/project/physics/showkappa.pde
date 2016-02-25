
# load geometry
geometry = d7_stokes.in2d                                        
# and mesh
mesh = d7_stokes.vol.gz

#load xfem-library and python-bindings
shared = libngsxfem_xfem                                       
shared = libngsxfem_xstokes                                       
# pymodule = d7_stokes

define constant heapsize = 1e9

define constant one = 1.0
define constant eps2 = -1e-8

define constant R = 0.05

define constant eps1 = 1e-1

define constant d = (0.2 + eps1)

define constant x0 = 0.0
define constant y0 = 0.0

define coefficient lset
(
 (        
  ((x-x0+d) < 0) 
   *
   (        
    ((y-y0+d) < 0) 
     * (sqrt((x-x0+d)*(x-x0+d)+(y-y0+d)*(y-y0+d))-R)
    +
    ((y-y0+d) > 0) * ((y-y0-d) < 0) 
     * (x0-x-d-R) 
    +
    ((y-y0-d) > 0) 
     * (sqrt((x-x0+d)*(x-x0+d)+(y-y0-d)*(y-y0-d))-R)
   )
  +
  ((x-x0+d) > 0) * ((x-x0-d) < 0) 
   *
   (        
    ((y-y0+d) < 0) 
     * (y0-y-d-R) 
    +
    ((y-y0+d) > 0) * ((y-y0-d) < 0) 
     * (-R) 
    +
    ((y-y0-d) > 0) 
     * (y-y0-d-R) 
   )
  +
  ((x-x0-d) > 0) 
   *
   (        
    ((y-y0+d) < 0) 
     * (sqrt((x-x0-d)*(x-x0-d)+(y-y0+d)*(y-y0+d))-R)
    +
    ((y-y0+d) > 0) * ((y-y0-d) < 0) 
     * (x-x0-d-R) 
    +
    ((y-y0-d) > 0) 
     * (sqrt((x-x0-d)*(x-x0-d)+(y-y0-d)*(y-y0-d))-R)
   )
 )
)

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

define constant zero = 0.0
define constant one = 1.0
define constant mu1 = 1.0
define constant mu2 = 1.0
define constant lambda = 1000.0
define constant delta = 1.0

define coefficient s
0,1,0,0,

define coefficient ghost
-1.0 ,

define linearform f -fespace=fescomp
xsource one zero -comp=2
xLBmeancurv one # naiv Laplace-Beltrami discretization 

define bilinearform a -fespace=fescomp -symmetric -linearform=f -printelmat
xstokes mu1 mu2
myghostpenalty ghost -comp=3

define preconditioner c -type=direct -bilinearform=a -inverse=pardiso #-test 

numproc bvp npbvp -gridfunction=uvp -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=100 -prec=1e-6

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

numproc xtonegpos npxtonegposu -xstd_gridfunction=uvp.1 -negpos_gridfunction=gf_u_negpos
numproc xtonegpos npxtonegposv -xstd_gridfunction=uvp.2 -negpos_gridfunction=gf_v_negpos
numproc xtonegpos npxtonegposp -xstd_gridfunction=uvp.3 -negpos_gridfunction=gf_p_negpos

#numproc vtkoutput npout -filename=kappa_
#        -coefficients=[lset]
#        -gridfunctions=[gf_u_negpos.1,gf_u_negpos.2,gf_v_negpos.1,gf_v_negpos.2,gf_p_negpos.1,gf_p_negpos.2]
#        -fieldnames=[levelset,uneg,upos,vneg,vpos,pneg,ppos]
#        -subdivision=1
