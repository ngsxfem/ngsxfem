# load geometry
geometry = square.in2d
# and mesh
#mesh = square.vol.gz
#mesh = square_trigs_2x2.vol.gz
mesh = square_trigs.vol.gz
#mesh = square_quad_coarse.vol.gz

shared = libngsxfem_xfem

define constant heapsize = 1e7

define fespace fesh1
       -type=spacetimefes 
       -type_space=h1ho
       -order_space=1
       -order_time=1

define fespace fesx
       -type=xfespace
       -spacetime
       -t0=0.0
       -t1=1.0
       -levelset=(x-y+z-0.45)

numproc informxfem npix 
        -fespace=fesh1
        -xfespace=fesx

define fespace fescomp
       -type=compound
       -spaces=[fesh1,fesx]

define gridfunction u -fespace=fescomp

define coefficient lset_past
(x-y-0.45),       

define coefficient lset_future
(x-y+1-0.45),       
       
#numproc shapetester npst -gridfunction=u

define bilinearform evalx_past -fespace=fescomp -nonassemble
xvis_st_past lset_past

define bilinearform evalx_future -fespace=fescomp -nonassemble
xvis_st_future lset_future

numproc drawflux npdf_past -solution=u -bilinearform=evalx_past -label=u_past -applyd
numproc drawflux npdf_future -solution=u -bilinearform=evalx_future -label=u_future -applyd

define constant one = 1.0
define constant zero = 0.0

define bilinearform a -fespace=fescomp -printelmat -print
stxmass one one

define linearform f -fespace=fescomp -print
stxsource one zero

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=direct -print


numproc visualization npviz -scalarfunction=u_future -comp=0

