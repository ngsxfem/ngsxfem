# load geometry
geometry = square.in2d
# and mesh
#mesh = square.vol.gz
#mesh = square_trigs_2x2.vol.gz
mesh = square_trigs.vol.gz
#mesh = square_quad_coarse.vol.gz

shared = libngsxfem_xstokes

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

define constant x0 = 0.5
define constant y0 = 0.5

define constant R = 0.33333333

define coefficient bconvneg
(bneg*wx,bpos*wy),

define coefficient bconvpos
(bneg*wx,bpos*wy),

define fespace fesu
       -type=spacetimefes 
       -type_space=h1ho
       -order_space=2
       -order_time=0
#       -dirichlet=[1,2]

define fespace fesp
       -type=spacetimefes 
       -type_space=h1ho
       -order_space=1
       -order_time=0
#       -dirichlet=[1,2]

define coefficient lset
((x-wx*z-x0)*(x-wx*z-x0)+(y-wy*z-y0)*(y-wy*z-y0)-R*R),
#(z),

define fespace fespx
       -type=xfespace
       -spacetime
       -t0=0.0
       -t1=0.05
       # -levelset=(x-y+z-0.375)
       # -levelset=((x-1.0*z-0.25)*(x-1.0*z-0.25)+(y-1.0*z-0.25)*(y-1.0*z-0.25)-0.04)
       -vmax=1.0
       -ref_space=0
       -ref_time=0

define fespace fescl 
       -type=lsetcontfespace

define fespace fespnegpos
       -type=compound
       -spaces=[fesp,fesp,fescl]

numproc informxfem npix 
        -fespace=fesp
        -xfespace=fespx
        -lsetcontfespace=fescl
        -coef_levelset=lset

define fespace fescomp_p
       -type=compound
       -spaces=[fesp,fespx]

define fespace fescomp
       -type=compound
       -spaces=[fesu,fesu,fescomp_p]
#       -dgjumps

define gridfunction u -fespace=fescomp
define gridfunction p_vis -fespace=fespnegpos

#numproc shapetester npst -gridfunction=u

define bilinearform evalu_past -fespace=fescomp -nonassemble
STtracepast one -comp=1
define bilinearform evalv_past -fespace=fescomp -nonassemble
STtracepast one -comp=2
define bilinearform evalu_future -fespace=fescomp -nonassemble
STtracefuture one -comp=1
define bilinearform evalv_future -fespace=fescomp -nonassemble
STtracefuture one -comp=2
define bilinearform evalpx_future -fespace=fescomp -nonassemble
stxvis_future one -comp=3
define bilinearform eval_pnegpos -fespace=fespnegpos -nonassemble
st_np_vis_future one
define bilinearform eval_pnegpos_past -fespace=fespnegpos -nonassemble
st_np_vis_past one


numproc drawflux npdfu_past -solution=u -bilinearform=evalu_past -label=u_past -applyd
numproc drawflux npdfv_past -solution=u -bilinearform=evalv_past -label=v_past -applyd
numproc drawflux npdfu_future -solution=u -bilinearform=evalu_future -label=u_future -applyd
numproc drawflux npdfv_future -solution=u -bilinearform=evalv_future -label=v_future -applyd
numproc drawflux npdfp_np -solution=p_vis -bilinearform=eval_pnegpos -label=p_negpos -applyd
numproc drawflux npdfp_np_test -solution=u -bilinearform=evalpx_future -label=p_test -applyd


# define coefficient rhsneg
# ((bneg)*1),

# define coefficient rhspos
# ((bpos)*1),

# define constant delta = 1

define bilinearform a -fespace=fescomp # -printelmat -print
STmass one told tnew -comp=1
STmass one told tnew -comp=2
# STtimeder one -comp=1
# STtimeder one -comp=2
stx_mass bneg bpos told tnew -comp=3
stx_stokes one one told tnew
#stx_laplace abneg abpos told tnew
# stx_nitsche_hansbo aneg apos bneg bpos lambda told tnew
# stx_timeder bneg bpos
# stx_convection bconvneg bconvpos told tnew
# stx_tracemass_past bneg bpos
# stx_lo_ghostpenalty abneg abpos told tnew delta
#stx_robin bneg_pen bpos_pen told tnew



define linearform f -fespace=fescomp # -print
STsource one told tnew -comp=1
STsource one told tnew -comp=2
stx_source one one told tnew -comp=3
#stx_tracesource_past binineg binipos
#stx_neumann bneg_bndneg_pen bpos_bndpos_pen told tnew

# numproc setvaluesx npsvx -gridfunction=u 
#         -coefficient_neg=bndneg 
#         -coefficient_pos=bndpos 
#         -told=0.0
#         -tnew=0.05
#         -boundary -print

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=direct # -print


# define coefficient veczero
# (0,0),

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

# numproc shapetester npst -gridfunction=u
