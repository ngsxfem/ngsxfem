
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
       -type=h1ho
       -order=2
#       -dirichlet=[1,2,3,4]

define fespace fesx
       -type=xfespace
       -levelset=((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)-0.09)

define coefficient lset
((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)-0.09),       
       
numproc informxfem npix 
        -fespace=fesh1
        -xfespace=fesx

define fespace fescomp
       -type=compound
       -spaces=[fesh1,fesx]

define gridfunction u -fespace=fescomp

numproc shapetester npst -gridfunction=u

define bilinearform evalx -fespace=fescomp -nonassemble
visx lset

numproc drawflux npdf -solution=u -bilinearform=evalx -label=utry -applyd

numproc visualization npviz -scalarfunction=utry -comp=0
