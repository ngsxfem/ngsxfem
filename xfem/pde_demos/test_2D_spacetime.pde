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
define constant bpos = 0.8

define constant aneg = 1.0
define constant apos = 1.0

define constant abneg = (aneg*bneg)
define constant abpos = (apos*bpos)

define constant lambda = 25.0

define constant told = 0.0
define constant tnew = 0.05

define constant wx = 0.5
define constant wy = 0.5

define constant binineg = 1.0
define constant binipos = 0.0

define constant pen = 1e7

define constant bneg_pen = (bneg*pen)
define constant bpos_pen = (bpos*pen)

define constant bndneg = (0.8*0.0)
define constant bndpos = (1.0*0.0)

define constant bneg_bndneg_pen = (bneg*pen*bndneg)
define constant bpos_bndpos_pen = (bpos*pen*bndpos)

define constant x0 = 0.5 #445
define constant y0 = 0.5

define constant R = 0.33333333

define coefficient bconvneg
(bneg*wx,bpos*wy),

define coefficient bconvpos
(bneg*wx,bpos*wy),

define coefficient lset
((x-wx*z-x0)*(x-wx*z-x0)+(y-wy*z-y0)*(y-wy*z-y0)-R*R),

define fespace fescomp
       -type=xh1fespace
       -type_space=h1ho
       -order_space=1
       -order_time=1
       -dirichlet=[1,2]
       -spacetime
       -t0=0.0
       -t1=0.05
       -vmax=1.0
       -ref_space=0
       -ref_time=0
#       -dgjumps

numproc informxfem npix 
        -xh1fespace=fescomp
        -coef_levelset=lset


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

define constant delta = 1

define bilinearform a -fespace=fescomp # -printelmat -print
#stx_mass bneg bpos told tnew
stx_laplace abneg abpos told tnew
stx_nitsche_hansbo aneg apos bneg bpos lambda told tnew
stx_timeder bneg bpos
stx_convection bconvneg bconvpos told tnew
stx_tracemass_past bneg bpos
#stx_lo_ghostpenalty abneg abpos told tnew delta
#stx_robin bneg_pen bpos_pen told tnew



define linearform f -fespace=fescomp # -print
#stx_source rhsneg rhspos told tnew
stx_tracesource_past binineg binipos
#stx_neumann bneg_bndneg_pen bpos_bndpos_pen told tnew

#numproc setvalues npsv -gridfunction=u.1 -coefficient=one 
#-boundary

numproc setvaluesx npsvx -gridfunction=u 
        -coefficient_neg=bndneg 
        -coefficient_pos=bndpos 
        -told=0.0
        -tnew=0.05
        -boundary -print

define preconditioner c -type=local -bilinearform=a -test -block
#define preconditioner c -type=bddc -bilinearform=a -test -block

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=gmres -preconditioner=c -maxsteps=1000 -prec=1e-6 # -print

# define coefficient veczero
# (0,0),

numproc calccond npcc -bilinearform=a -inverse=pardiso #-symmetric

# numproc xdifference npxd -solution=u 
#         -solution_n=two
#         -solution_p=one
#         -derivative_n=veczero
#         -derivative_p=veczero
#         -levelset=lset
#         -interorder=2
#         -henryweight_n=1.0
#         -henryweight_p=1.0


numproc visualization npviz -scalarfunction=u_future -subdivision=3

