
# load geometry
geometry = d1_approx.in2d                                        
# and mesh
#mesh = d1_approx.vol.gz
mesh = d1_approx2.vol.gz

#load xfem-library and python-bindings
shared = libngsxfem_xfem                                       
shared = libngsxfem_py                                       

define constant heapsize = 1e9
constant R = 0
        
# interface description as zero-level
define coefficient lset
( sqrt(x*x+y*y) - R),

fespace fes_p1 -type=h1ho -order=1
gridfunction lset_p1 -fespace=fes_p1
numproc setvalues npsv -gridfunction=lset_p1 -coefficient=lset
        
# evaluate l2-error (difference between prescribed solution and numerical solution)

fespace fes_deform -type=h1ho -order=2 -dim=2
gridfunction deform -fespace=fes_deform
                
numproc xgeomtest npxd 
        -levelset=lset
        -deformation=deform

numproc draw npdr -coefficient=lset -label=levelset



        
