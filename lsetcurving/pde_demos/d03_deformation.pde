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
                        
### project level set function into a finite element space of order k
fespace fes_ho -type=h1ho -order=3
gridfunction lset_ho -fespace=fes_ho

numproc setvalues npsv -gridfunction=lset_ho -coefficient=lset

### project this level set function into a finite element space of order 1
fespace fes_p1 -type=h1ho -order=1
gridfunction lset_p1 -fespace=fes_p1

numproc interpolatep1 npipp1b -gridfunction_ho=lset_ho -gridfunction_p1=lset_p1

### determine the deformation 
fespace fes_deform -type=h1ho -order=3 -vec
        
gridfunction deform -fespace=fes_deform

constant lset_lower_bound = 0.0
constant lset_upper_bound = 0.0
constant one = 1.0
constant threshold = 0.1
                                        
bilinearform mdeform -fespace=fes_deform -symmetric
restrictedmass_4 one lset_p1 lset_lower_bound lset_upper_bound --comp=1        
restrictedmass_4 one lset_p1 lset_lower_bound lset_upper_bound --comp=2
#mass 1.0 --comp=1        

linearform fdeform -fespace=fes_deform
shiftsource_5 lset_p1 lset_ho threshold lset_lower_bound lset_upper_bound

define preconditioner cdeform -type=direct -bilinearform=mdeform
              
numproc bvp npbvpdeform -gridfunction=deform -bilinearform=mdeform -linearform=fdeform -solver=cg -preconditioner=cdeform -maxsteps=1000 -prec=1e-6 # -print
                       

          
                