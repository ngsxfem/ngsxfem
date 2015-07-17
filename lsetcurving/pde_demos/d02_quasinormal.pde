geometry = d01_testgeom.in2d
mesh = d01_testgeom.vol.gz
# mesh = d99_testgeom_unstr.vol.gz
shared = libngsxfem_xfem
shared = libngsxfem_lsetcurving

#define constant geometryorder = 8
        
define constant heapsize = 1e9
constant R = 2.0/3.0

define coefficient lset
(sqrt((x)*(x)+(y)*(y)) - R),
                        
fespace fes_p1 -type=h1ho -order=1
gridfunction lset_p1 -fespace=fes_p1

fespace fes_ho -type=h1ho -order=3
gridfunction lset_ho -fespace=fes_ho

numproc setvalues npsv -gridfunction=lset_ho -coefficient=lset

fespace fes_normal -type=h1ho -order=3 -vec #-dirichlet=[1,2,3,4,5,6]
gridfunction qn -fespace=fes_normal
        
bilinearform mqn -fespace=fes_normal
mass 1.0 --comp=1        
mass 1.0 --comp=2
#mass 1.0 --comp=1        

coefficient test1
(grad_lset_ho*(1,0))

coefficient test2
(grad_lset_ho*(0,1))
                
linearform fqn -fespace=fes_normal
source test1 --comp=1
source test2 --comp=2
# source (grad_lset_ho*(0,1)) --comp=2


define preconditioner cqn -type=direct -bilinearform=mqn
              
numproc bvp npbvpqn -gridfunction=qn -bilinearform=mqn -linearform=fqn -solver=cg -preconditioner=cqn -maxsteps=1000 -prec=1e-6 # -print
                       

          
                