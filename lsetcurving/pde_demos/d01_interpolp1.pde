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
gridfunction lset_p1_a -fespace=fes_p1
gridfunction lset_p1_b -fespace=fes_p1

fespace fes_ho -type=h1ho -order=4
gridfunction lset_ho -fespace=fes_ho

numproc setvalues npsv -gridfunction=lset_ho -coefficient=lset

#interpolate the coefficient function        
numproc interpolatep1 npipp1a -coefficient=lset -gridfunction_p1=lset_p1_a
#interpolate the gridfunction lset_ho        
numproc interpolatep1 npipp1b -gridfunction_ho=lset_ho -gridfunction_p1=lset_p1_b


