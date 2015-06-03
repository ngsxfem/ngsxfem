
# load geometry
geometry = d1_approx.in2d                                        
# and mesh
mesh = d1_approx.vol.gz

#load xfem-library and python-bindings
shared = libngsxfem_xfem                                       
shared = libngsxfem_py                                       

define constant heapsize = 1e9

define constant R = 0.4
define constant one = 1.0

# interface description as zero-level
define coefficient lset
( sqrt(x*x+y*y) - R),

# solution u- in 'inner' domain (Omega_1)
define coefficient solneg
0.5,

# solution u+ in 'outer' domain (Omega_2)
define coefficient solpos
(sin((x*x+y*y-R*R))),

# use an "extended" continuous finite element space
# you may change the order here
define fespace fescomp                                            
       -type=xstdfespace                                      
       -type_std=h1ho                                          
       -order=1                                                  
       -ref_space=1                                           
#       -empty
               
        
#update "extended" part of XFE space:
numproc informxfem npix                                     
        -xstdfespace=fescomp                                      
        -coef_levelset=lset

define gridfunction u -fespace=fescomp

# integration on sub domains
define linearform f -fespace=fescomp
xsource solneg solpos                             

# integration on sub domains
define bilinearform a -fespace=fescomp -symmetric -linearform=f
xmass 1.0 1.0

#define preconditioner c -type=local -bilinearform=a -test #-block           
define preconditioner c -type=direct -bilinearform=a -inverse=pardiso #-test 

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6

numproc visualization npviz -scalarfunction=u 
    -minval=0 -maxval=1 
    -nolineartexture -deformationscale=1 -subdivision=4

# evaluate l2-error (difference between prescribed solution and numerical solution)
numproc xdifference npxd 
        -solution=u 
        -solution_n=solneg
        -solution_p=solpos
        -intorder=3



define fespace fesstd
        -type=h1ho
        -order=1

define fespace fesnegpos
        -type=compound
        -spaces=[fesstd,fesstd]

define gridfunction gf_negpos -fespace=fesnegpos

numproc xtonegpos npxtonegpos -xstd_gridfunction=u -negpos_gridfunction=gf_negpos

numproc vtkoutput npout -filename=d5_output
        -coefficients=[lset]
        -gridfunctions=[gf_negpos.1,gf_negpos.2]
        -fieldnames=[levelset,uneg,upos]
        -subdivision=1
