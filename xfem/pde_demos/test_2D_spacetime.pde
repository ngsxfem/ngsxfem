# load geometry
geometry = square.in2d
# and mesh
#mesh = square.vol.gz
#mesh = square_trigs_2x2.vol.gz
mesh = square_trigs.vol.gz
#mesh = square_quad_coarse.vol.gz

shared = libngsxfem_xfem

define constant heapsize = 1e7

define constant zero = 0.0
define constant one = 1.0

define constant two = 2.0

define constant bneg = 2.0
define constant bpos = 1.0

define constant aneg = 10.0
define constant apos = 10.0

define constant abneg = (aneg*bneg)
define constant abpos = (apos*bpos)

define constant lambda = 10.0

define constant told = 0.0
define constant tnew = 0.005

define constant wx = 1.0
define constant wy = 1.0

define constant binineg = 1.0
define constant binipos = 0.0

define coefficient bconvneg
(bneg*wx,bpos*wy),

define coefficient bconvpos
(bneg*wx,bpos*wy),

define fespace fesh1
       -type=spacetimefes 
       -type_space=h1ho
       -order_space=1
       -order_time=1
#       -dirichlet=[1,2]

define coefficient lset
((x-1.0*z-0.25)*(x-1.0*z-0.25)+(y-1.0*z-0.25)*(y-1.0*z-0.25)-0.04),
#(z),

define fespace fesx
       -type=xfespace
       -spacetime
       -t0=0.0
       -t1=0.005
       # -levelset=(x-y+z-0.375)
       -levelset=((x-1.0*z-0.25)*(x-1.0*z-0.25)+(y-1.0*z-0.25)*(y-1.0*z-0.25)-0.04)

numproc informxfem npix 
        -fespace=fesh1
        -xfespace=fesx
        -coef_levelset=lset

define fespace fescomp
       -type=compound
       -spaces=[fesh1,fesx]

define gridfunction u -fespace=fescomp
       
#numproc shapetester npst -gridfunction=u

define bilinearform evalx_past -fespace=fescomp -nonassemble
stxvis_past one

define bilinearform evalx_future -fespace=fescomp -nonassemble
stxvis_future one

numproc drawflux npdf_past -solution=u -bilinearform=evalx_past -label=u_past -applyd
numproc drawflux npdf_future -solution=u -bilinearform=evalx_future -label=u_future -applyd

define coefficient rhsneg
((bneg)*1),
#(sin(x)),

define coefficient rhspos
((bpos)*1),
#(cos(x)),

define bilinearform a -fespace=fescomp # -printelmat -print
#stx_mass bneg bpos told tnew
stx_laplace abneg abpos told tnew
stx_nitsche_halfhalf aneg apos bneg bpos lambda told tnew
stx_timeder bneg bpos
stx_convection bconvneg bconvpos told tnew
stx_tracemass_past bneg bpos



define linearform f -fespace=fescomp # -print
#stx_source rhsneg rhspos told tnew
stx_tracesource_past binineg binipos


#numproc setvalues npsv -gridfunction=u.1 -coefficient=one 
#-boundary


numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=direct # -print


numproc visualization npviz -scalarfunction=u_future -subdivision=3

