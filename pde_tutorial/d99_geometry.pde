
geometry = d1_approx.in2d                                        
mesh = d1_approx.vol.gz

# # load geometry
# geometry = d1_approx.in2d                                        
# # and mesh
# mesh = d1_approx.vol.gz
# #mesh = d1_approx2.vol.gz

#load xfem-library and python-bindings
shared = libngsxfem_xfem                                       
shared = libngsxfem_py                                       

define constant heapsize = 1e9
constant R = 0.5
        
#interface description as zero-level

# define coefficient lset
# # # ( x*x+y*y - R*R),
# ( sqrt(x*x+y*y) - R),

# general constants
define constant zero = 0.0
define constant one = 1.0
define constant two = 2.0

define constant C1 = 1.0
define constant C2 = 4.0

#geometry constants
define constant r0 = 0.5+1e-12
define constant omega = 5
define constant omega2 = (1.2*pi) 

define coefficient lset
(
 sqrt(x*x+y*y)-(r0+0.2*sin(omega*atan2(x,y)))
)

# define coefficient lset
# ( sqrt(x*x+4*y*y) - R),
                
# define constant R = 0.05

# define constant d = 0.2

# define constant x0 = 0.5
# define constant y0 = 0.5
        
# define coefficient lset
# (
#  (        
#   ((x-x0+d) < 0) 
#    *
#    (        
#     ((y-y0+d) < 0) 
#      * (sqrt((x-x0+d)*(x-x0+d)+(y-y0+d)*(y-y0+d))-R)
#     +
#     ((y-y0+d) > 0) * ((y-y0-d) < 0) 
#      * (x0-x-d-R) 
#     +
#     ((y-y0-d) > 0) 
#      * (sqrt((x-x0+d)*(x-x0+d)+(y-y0-d)*(y-y0-d))-R)
#    )
#   +
#   ((x-x0+d) > 0) * ((x-x0-d) < 0) 
#    *
#    (        
#     ((y-y0+d) < 0) 
#      * (y0-y-d-R) 
#     +
#     ((y-y0+d) > 0) * ((y-y0-d) < 0) 
#      * (-R) 
#     +
#     ((y-y0-d) > 0) 
#      * (y-y0-d-R) 
#    )
#   +
#   ((x-x0-d) > 0) 
#    *
#    (        
#     ((y-y0+d) < 0) 
#      * (sqrt((x-x0-d)*(x-x0-d)+(y-y0+d)*(y-y0+d))-R)
#     +
#     ((y-y0+d) > 0) * ((y-y0-d) < 0) 
#      * (x-x0-d-R) 
#     +
#     ((y-y0-d) > 0) 
#      * (sqrt((x-x0-d)*(x-x0-d)+(y-y0-d)*(y-y0-d))-R)
#    )
#  )
# )


fespace fes_p1 -type=h1ho -order=1
gridfunction lset_p1 -fespace=fes_p1

fespace fes_p2 -type=h1ho -order=2
gridfunction lset_p2 -fespace=fes_p2
                
fespace fes_deform -type=h1ho -order=2 -dim=2
gridfunction deform -fespace=fes_deform
                
numproc xgeomtest npxd 
        -gf_levelset_p1=lset_p1
        -gf_levelset_p2=lset_p2
        -levelset=lset
        -gf_levelset=lset_p2
        -deformation=deform

numproc draw npdr -coefficient=lset -label=levelset


numproc visualization npvis -scalarfunction=lset_p1 -vectorfunction=deform -deformationscale=1 -subdivision=5 -minval=0 -maxval=0



        
