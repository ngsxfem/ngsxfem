
# load geometry
geometry = square.in2d

# and mesh
#mesh = square.vol.gz
#mesh = square_trigs_2x2.vol.gz
mesh = square_trigs.vol.gz
#mesh = square_quad_coarse.vol.gz

shared = libngsxfem_xfem

define constant heapsize = 1e8

define constant zero = 0.0
define constant one = 1.0

define constant two = 2.0

define constant pen = 1e7

define constant x0 = 0.5
define constant y0 = 0.5

define constant bneg = 1.0
define constant bpos = 2.0

define constant aneg = 0.1
define constant apos = 0.1

define constant abneg = (aneg*bneg)
define constant abpos = (apos*bpos)

define constant bneg_pen = (pen*bneg)
define constant bpos_pen = (pen*bpos)

define constant bneg_bndvalneg_pen = (pen*bneg*bpos)
define constant bpos_bndvalpos_pen = (pen*bpos*bneg)

define constant lambda = 15.0

define constant R = 0.1251 #0.33333333

define fespace fesh1
       -type=h1ho
       -order=1
       -dirichlet=[1,2]

define fespace fesx
       -type=xfespace
#       -levelset=(x-0.55)
#       -levelset=((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)-0.09)

define coefficient lset
#(x-0.55),
(abs(x-x0)-R),
#((x-x0)*(x-x0)+(y-y0)*(y-y0)-R*R),       
       
numproc informxfem npix 
        -fespace=fesh1
        -xfespace=fesx
        -coef_levelset=lset

define fespace fescomp
       -type=compound
       -spaces=[fesh1,fesx]
       -dgjumps

define gridfunction u -fespace=fescomp

#numproc shapetester npst -gridfunction=u

define bilinearform evalx -fespace=fescomp -nonassemble
xvis one

numproc drawflux npdf -solution=u -bilinearform=evalx -label=utry -applyd

define bilinearform a -fespace=fescomp # -printelmat -print
#xmass one one
xlaplace abneg abpos
xnitsche_halfhalf aneg apos bneg bpos lambda
#xrobin bneg_pen bpos_pen
lo_ghostpenalty aneg apos

define linearform f -fespace=fescomp # -print
xsource bneg bpos
#xneumann bneg_bndvalneg_pen bpos_bndvalpos_pen

#numproc setvalues npsv -gridfunction=u.1 -coefficient=zero -boundary

numproc setvaluesx npsvx -gridfunction=u -coefficient_neg=two -coefficient_pos=one -boundary -print

define preconditioner c -type=local -bilinearform=a
#define preconditioner c -type=direct -bilinearform=a

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=gmres -preconditioner=c # -print

define coefficient veczero
(0,0),

# numproc xdifference npxd -solution=u 
#         -function_n=two
#         -function_p=one
#         -derivative_n=veczero
#         -derivative_p=veczero
#         -levelset=lset
#         -interorder=2
#         -henryweight_n=1.0
#         -henryweight_p=1.0

numproc visualization npviz -scalarfunction=utry #-comp=0
