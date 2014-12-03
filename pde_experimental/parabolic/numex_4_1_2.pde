#
# solve the Poisson equation -Delta u = f
#
# with boundary conditions
#      u = 0  on Gamma1
#  du/dn = 1  on Gamma2

# load geometry
geometry = numex_4_1_1.in2d

# and mesh
mesh = numex_4_1_1_reg.vol.gz

shared = libngsxfem_spacetime
shared = libngsxfem_xfem
shared = libngsxfem_parabolic

define constant heapsize = 1e8

define constant one = 1.0

define constant aneg = 1.0
define constant apos = 2.0

define constant bneg = 1.5
define constant bpos = 1.0

define constant wx = 1.0
define constant wy = 0.0

define constant R = (1.0/3.0)

# define constant b = (0.5/(R*R)*(apos/aneg*pi*cos(pi*R)-bpos/bneg*sin(pi*R)/R))
# define constant a = (bpos/bneg*sin(pi*R)/R - b*R*R)

define constant b = 6.34294
define constant a = 1.02728

define constant k = 1.5

define coefficient solneg
(sin(k*pi*z)*( a*(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))
              +b*(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))
                *(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))
                *(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y))))),


define coefficient solpos
(sin(k*pi*z) * sin(pi*(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y))))),

define coefficient bconvneg
((bneg*0.5*cos(2.0*pi*z)), 0.0),

define coefficient bconvpos
((bpos*0.5*cos(2.0*pi*z)), 0.0),

define coefficient binineg
(bneg*sin(k*pi*z)*( a*(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))
                   +b*(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))
                     *(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))
                     *(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y))))),


define coefficient binipos
(bpos*sin(k*pi*z) * sin(pi*(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y))))),

define coefficient brhsneg
(bneg*
      (k*pi*cos(k*pi*z)*(a*(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))
                        +b*(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))
                          *(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))
                          *(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y))))
       -aneg*sin(k*pi*z)*( (1+y*y*(1-y)*(1-y)*(2-y)*(2-y))*(6*b*(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y))))
                          -((2-6*y+3*y*y)*(a+3*b*(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))
                                                *(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))
                                          )
                           )
                         )
      )
),

define coefficient brhspos
(bpos*
      (k*pi*cos(k*pi*z)*(sin(pi*(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))))
       -apos*sin(k*pi*z)*( (1+y*y*(1-y)*(1-y)*(2-y)*(2-y))*(-pi*pi*sin(pi*(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))))
                          -((2-6*y+3*y*y)*(pi*cos(pi*(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y))))))
                         )
      )
),

define fespace fesh1
       -type=spacetimefes 
       -type_space=h1ho
       -order_space=1
       -order_time=1
       -dirichlet=[2,4]

define coefficient coef_lset
(abs(x-(0.25/pi*sin(2*pi*z))-(7.0/8.0+0.25*y*y*(2-y)*(2-y)))-R),

define fespace fesx
       -type=xfespace
       -spacetime
       -vmax=0.25
       -ref_space=0
       -ref_time=0

define fespace fescl 
       -type=lsetcontfespace

define fespace fesnegpos
       -type=compound
       -spaces=[fesh1,fesh1,fescl]

numproc informxfem npix 
        -fespace=fesh1
        -xfespace=fesx
        -lsetcontfespace=fescl
        -coef_levelset=coef_lset

define fespace fescomp
       -type=compound
       -spaces=[fesh1,fesx]
#       -dgjumps


define gridfunction u -fespace=fescomp
define gridfunction u_vis -fespace=fesnegpos

numproc stx_solveinstat npsi 
        -initialneg=binineg
        -initialpos=binipos
        -beta_conv_neg=bconvneg
        -beta_conv_pos=bconvpos
        -beta_rhs_neg=brhsneg
        -beta_rhs_pos=brhspos
        -beta_ini_neg=binineg
        -beta_ini_pos=binipos
        -boundary_neg=solneg
        -boundary_pos=solpos
        -gf_vis=u_vis
        -gridfunction=u
        -solver=pardiso 
        -fespace=fescomp
        -fespacevis=fesnegpos
        -dt=0.1
        -tstart=0
        -tend=1.0
        -aneg=1
        -apos=2
        -bneg=1.5
        -bpos=1.0
        -lambda=20.0
#        -pause_after_step=2000
        -solution_n=solneg
        -solution_p=solpos
        -levelset=coef_lset
#        -userstepping
#        -calccond
#        -ghostpenalty
#        -delta=0.5
#        -direct
#         -minimal_stabilization

define coefficient veczero
(0,0),

numproc xdifference npxd -solution=u 
        -solution_n=solneg
        -solution_p=solpos
        # -derivative_n=veczero
        # -derivative_p=veczero
        -levelset=coef_lset
        -intorder=4
        -henryweight_n=1.5
        -henryweight_p=1.0
        -time=1.0 # coeff function not yet of higher dimension...



# define bilinearform evalx_past -fespace=fescomp -nonassemble
# stxvis_past one

define bilinearform evalx_future -fespace=fescomp -nonassemble
stxvis_future one

# define bilinearform evalx_neg -fespace=fesnegpos -nonassemble
# STtracefuture one -comp=1

# define bilinearform evalx_pos -fespace=fesnegpos -nonassemble
# STtracefuture one -comp=2

# numproc drawflux npdf_past -solution=u -bilinearform=evalx_past -label=u_past -applyd
numproc drawflux npdf_future -solution=u -bilinearform=evalx_future -label=u_future -applyd

# numproc drawflux npdf_past -solution=u_vis -bilinearform=evalx_neg -label=u_neg -applyd
# numproc drawflux npdf_future -solution=u_vis -bilinearform=evalx_pos -label=u_pos -applyd

# numproc visualization npviz 
#         -scalarfunction=u_future 
#         -subdivision=3
#         -nolineartexture        
        # -deformationscale=0.3 
        # -subdivision=2 
        # -minval=0.0 
        # -maxval=1.0

define bilinearform eval_negpos -fespace=fesnegpos -nonassemble
st_np_vis_future one

define bilinearform eval_negpos_past -fespace=fesnegpos -nonassemble
st_np_vis_past one

numproc drawflux npdf_np -solution=u_vis -bilinearform=eval_negpos -label=u_negpos -applyd
numproc drawflux npdf_np_past -solution=u_vis -bilinearform=eval_negpos_past -label=u_negpos_past -applyd

numproc visualization npvis -scalarfunction=u_negpos -nolineartexture -deformationscale=1.0 -subdivision=0 
        -minval=-1.0 -maxval=1.0
