########################################################################################
# solve poisson equation
#    - eps Delta u = f in Omega,  u=g on Gamma0
########################################################################################


geometry = square.in2d
mesh = square.vol.gz

# geometry = cube.geo
# mesh = cube.vol

shared = libngsxfem_meancurv
define constant heapsize = 10000000


# ho-fespace: compound space for (u, lambda) ######
define fespace v -type=h1ho -order=5

define gridfunction u -fespace=v 

define constant one  = 1
define constant zero = 0

## boundary terms ########################################
# define coefficient initialvals
# (sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))-0.25),

# define coefficient initialvals
# (sqrt(0.6*(sin(x)-0.5)*(sin(x)-0.5)+(y-0.5)*(y-0.5))-0.25),
define coefficient initialvals
(
(sqrt((x-0.3)*(x-0.3)+(y-0.5)*(y-0.5))-0.25 < sqrt((x-0.7)*(x-0.7)+(y-0.5)*(y-0.5))-0.25) * (sqrt((x-0.3)*(x-0.3)+(y-0.5)*(y-0.5))-0.25) 
+
(sqrt((x-0.3)*(x-0.3)+(y-0.5)*(y-0.5))-0.25 > sqrt((x-0.7)*(x-0.7)+(y-0.5)*(y-0.5))-0.25) * (sqrt((x-0.7)*(x-0.7)+(y-0.5)*(y-0.5))-0.25)
),


define linearform f -fespace=v

define bilinearform a -fespace=v
meancurv_stiff

define bilinearform m -fespace=v
meancurv_mass

numproc setvalues npsv -gridfunction=u -coefficient=initialvals
#numproc setvalues npsv -gridfunction=ucip -coefficient=cdir -boundary

numproc meancurv np1 -bilinearforma=a -bilinearformm=m -gridfunction=u -dt=0.0001 -tend=0.07

numproc visualization vis1 -scalarfunction=u -subdivision=2 -minval=-1e-8 -maxval=1e-8 #-deformationscale=1


