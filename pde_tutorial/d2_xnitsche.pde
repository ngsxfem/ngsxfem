
# load geometry
geometry = d2_xnitsche.in2d

# and mesh
mesh = d2_xnitsche.vol.gz

#load xfem-library
shared = libngsxfem_xfem                                    
shared = libngsxfem_py                                      
pymodule = d1_approx

define constant heapsize = 1e8

# center of circle (domain 1)
define constant x0 = 0.1
define constant y0 = 0

# radius of circle (domain 1)
define constant R = 0.5

define coefficient lset
( sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)) - R),

numproc draw npd -coefficient=lset -label=levelset

# henry weights
define constant bneg = 3.0
define constant bpos = 4.0

define constant aneg = 0.3
define constant apos = 0.4

define coefficient rhspos
1,

define coefficient rhsneg
1,

define coefficient bndpos
0,

define coefficient bndneg
0,

define constant lambda = 2.0

define constant one = 1.0

define fespace fescomp                                            
       -type=xstdfespace                                      
       -type_std=h1ho                                          
       -order=1                                                  
       -dirichlet=[1,2]                                         
       -ref_space=1                                           
#       -empty                                                      
#       -dgjumps                                                  

numproc informxfem npix
        -xstdfespace=fescomp
        -coef_levelset=lset

define gridfunction u -fespace=fescomp

define linearform f -fespace=fescomp
xsource rhsneg rhspos

define bilinearform a -fespace=fescomp -symmetric -linearform=f
xlaplace aneg*bneg apos*bpos                                    
xnitsche_hansbo aneg apos bneg bpos lambda                     
# xnitsche_minstab_hansbo aneg apos bneg bpos           
#lo_ghostpenalty aneg apos one

numproc setvaluesx npsvx -gridfunction=u -coefficient_neg=bndneg -coefficient_pos=bndpos -boundary

#define preconditioner c -type=local -bilinearform=a -test #-block
define preconditioner c -type=direct -bilinearform=a -inverse=pardiso #-test

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6 

numproc visualization npviz -scalarfunction=u 
    -minval=0 -maxval=0.3
    -nolineartexture -deformationscale=1 -subdivision=4
