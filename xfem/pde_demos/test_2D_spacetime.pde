# load geometry
geometry = square.in2d
# and mesh
#mesh = square.vol.gz
#mesh = square_trigs_2x2.vol.gz
mesh = square_trigs.vol.gz
#mesh = square_quad_coarse.vol.gz

shared = libngsxfem_xfem

define constant heapsize = 1e7

define constant told = 0.0
define constant tnew = 0.1

define constant one = 1.0

define fespace fesh1
       -type=spacetimefes 
       -type_space=h1ho
       -order_space=1
       -order_time=1
#       -dirichlet=[1,2]

define fespace fesx
       -type=xfespace
       -spacetime
       -t0=0.0
       -t1=0.1
       # -levelset=(x-y+z-0.375)
       -levelset=((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)-0.09-z)

numproc informxfem npix 
        -fespace=fesh1
        -xfespace=fesx

define fespace fescomp
       -type=compound
       -spaces=[fesh1,fesx]

define gridfunction u -fespace=fescomp

define coefficient lset_past
#(x-y+told-0.375),       
((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)-0.09-told)

define coefficient lset_future
#(x-y+tnew-0.375),       
((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)-0.09-tnew)
       
#numproc shapetester npst -gridfunction=u

define bilinearform evalx_past -fespace=fescomp -nonassemble
stxvis_past one

define bilinearform evalx_future -fespace=fescomp -nonassemble
stxvis_future one

numproc drawflux npdf_past -solution=u -bilinearform=evalx_past -label=u_past -applyd
numproc drawflux npdf_future -solution=u -bilinearform=evalx_future -label=u_future -applyd

define coefficient rhs1
0,
#(sin(x)),

define coefficient rhs2
1,
#(cos(x)),

define bilinearform a -fespace=fescomp # -printelmat -print
stxmass one one told tnew

define linearform f -fespace=fescomp # -print
stxsource rhs1 rhs2 told tnew


#numproc setvalues npsv -gridfunction=u.1 -coefficient=one 
#-boundary


numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=direct # -print


numproc visualization npviz -scalarfunction=u_future

