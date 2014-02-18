
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

define constant bneg = 1.0
define constant bpos = 1.0

define constant aneg = 1.0
define constant apos = 50.0

define constant lambda = 500.0

define fespace fesh1
       -type=h1ho
       -order=1
       -dirichlet=[1,2]

define fespace fesx
       -type=xfespace
#       -levelset=(x-0.55)
       -levelset=((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)-0.09)

define coefficient lset
#(x-0.55),
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
xvis one

numproc drawflux npdf -solution=u -bilinearform=evalx -label=utry -applyd

define bilinearform a -fespace=fescomp # -printelmat -print
#xmass one one
xlaplace aneg apos
xnitsche_halfhalf aneg apos bneg bpos lambda

define linearform f -fespace=fescomp # -print
xsource one one

numproc setvalues npsv -gridfunction=u.1 -coefficient=zero -boundary

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=direct # -print

numproc visualization npviz -scalarfunction=utry #-comp=0
