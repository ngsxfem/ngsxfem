#
# solve the Poisson equation -Delta u = f
#
# with boundary conditions
#      u = 0  on Gamma1
#  du/dn = 1  on Gamma2

# load geometry
geometry = square.in2d

# and mesh
#mesh = square.vol.gz
mesh = square_trigs.vol.gz
#mesh = square_quad_coarse.vol.gz

#shared = libngsxfem_test
#shared = libngsxfem_common
shared = libngsxfem_parabolic

define constant heapsize = 1e7

define constant one = 1.0

define constant bneg = 1.0
define constant bpos = 2.0

define constant wx = 1.0
define constant wy = 1.0

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

define fespace fesx
       -type=xfespace
       -spacetime
       -t0=0.0
       -t1=0.0005
       # -levelset=(x-y+z-0.375)
       -levelset=((x-1.0*z-0.4)*(x-1.0*z-0.4)+(y-1.0*z-0.4)*(y-1.0*z-0.4)-0.04)

numproc informxfem npix 
        -fespace=fesh1
        -xfespace=fesx

define fespace fescomp
       -type=compound
       -spaces=[fesh1,fesx]

define gridfunction u -fespace=fescomp

numproc stx_solveinstat npsi 
        -initialneg=binineg
        -initialpos=binipos
        -beta_conv_neg=bconvneg
        -beta_conv_pos=bconvpos
        # -linearform=f 
        -gridfunction=u
        -solver=pardiso 
        -fespace=fescomp
        -dt=0.0005
        -tend=0.1
        -userstepping
        -aneg=10.0
        -apos=10.0
        -bneg=1.0
        -bpos=2.0
        -lambda=100.0

define bilinearform evalx_past -fespace=fescomp -nonassemble
stxvis_past one

define bilinearform evalx_future -fespace=fescomp -nonassemble
stxvis_future one

numproc drawflux npdf_past -solution=u -bilinearform=evalx_past -label=u_past -applyd
numproc drawflux npdf_future -solution=u -bilinearform=evalx_future -label=u_future -applyd

numproc visualization npviz 
        -scalarfunction=u_future 
        -subdivision=3
        -nolineartexture        
        # -deformationscale=0.3 
        # -subdivision=2 
        # -minval=0.0 
        # -maxval=1.0


numproc visualization npvis -scalarfunction=u_future -nolineartexture -deformationscale=0.3 -subdivision=2 #-minval=0.0 -maxval=1.0
