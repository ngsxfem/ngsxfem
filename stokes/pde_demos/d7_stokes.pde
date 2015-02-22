
# load geometry
geometry = d7_stokes.in2d                                        
# and mesh
mesh = d7_stokes.vol.gz

#load xfem-library and python-bindings
shared = libngsxfem_xfem                                       
shared = libngsxfem_xstokes                                       
# pymodule = d7_stokes

define constant heapsize = 1e9

define constant R = 0.4
define constant one = 1.0

# interface description as zero-level
define coefficient lset
( sqrt(x*x+y*y) - R),

# define coefficient lset
# ( x-0.127378 ),



define fespace fesuv
       -type=h1ho
       -order=2
       -dirichlet=[1,2,3,4]

define fespace fesuvx
       -type=xfespace
       -ref_space=0
       -empty

define fespace fesp
       -type=h1ho
       -order=1
       # -dirichlet=[1,2,3,4]

define fespace fespx
       -type=xfespace
       -ref_space=0
       # -empty

define fespace fescl 
       -type=lsetcontfespace

numproc informxfem npi_uvx 
        -fespace=fesuv
        -xfespace=fesuvx
        -lsetcontfespace=fescl
        -coef_levelset=lset

numproc informxfem npi_px 
        -fespace=fesp
        -xfespace=fespx
        -lsetcontfespace=fescl
        -coef_levelset=lset

define fespace fescomp_u
       -type=compound
       -spaces=[fesuv,fesuvx]

define fespace fescomp_v
       -type=compound
       -spaces=[fesuv,fesuvx]

define fespace fescomp_p
       -type=compound
       -spaces=[fesp,fespx]

define fespace fescomp
       -type=compound
       -spaces=[fescomp_u,fescomp_v,fescomp_p]

define gridfunction uvp -fespace=fescomp

define constant one = 1.0
define constant lambda = 10.0

define coefficient s
0,1,0,0,

numproc setvaluesx npsvx -gridfunction=uvp.2 -coefficient_neg=s -coefficient_pos=s -boundary

# integration on sub domains
define linearform f -fespace=fescomp
# xsource s s -comp=2

# integration on sub domains
define bilinearform a -fespace=fescomp -symmetric -linearform=f -printelmat
xstokes one one 
# xlaplace one one -comp=1
xnitsche one one one one lambda -comp=1
xnitsche one one one one lambda -comp=2
# xmass one one -comp=1
# xmass 1.0 1.0 -comp=2

#define preconditioner c -type=local -bilinearform=a -test #-block           
define preconditioner c -type=direct -bilinearform=a -inverse=pardiso #-test 

numproc bvp npbvp -gridfunction=uvp -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6

numproc visualization npviz -scalarfunction=u 
    -minval=0 -maxval=1 
    -nolineartexture -deformationscale=1 -subdivision=4
