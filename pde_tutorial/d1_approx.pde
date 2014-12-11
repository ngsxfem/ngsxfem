
# load geometry
geometry = d1_approx.in2d

# and mesh
mesh = d1_approx.vol.gz

#load xfem-library
shared = libngsxfem_xfem
shared = libngsxfem_py

pymodule = d1_approx

define constant heapsize = 1e9

define constant R = 0.4
define constant one = 1.0

# interface description as zero-level
define coefficient lset
( sqrt(x*x+y*y) - R),       

# solution u- in 'inner' domain (Omega_1)
define coefficient solneg
0.5,

# solution u+ in 'outer' domain (Omega_2)
define coefficient solpos
(sin((x*x+y*y-R*R))),

# use an "extended" continuous finite element space
# you may change the order here
define fespace fescomp
       -type=xstdfespace
       -type_std=h1ho
       -order=1
#       -ref_space=5
#       -empty

#update "extended" part of XFE space:
numproc informxfem npix 
        -xstdfespace=fescomp
        -coef_levelset=lset

define gridfunction u -fespace=fescomp

# integration on sub domains
define linearform f -fespace=fescomp
xsource solneg solpos

# integration on sub domains
define bilinearform a -fespace=fescomp -symmetric -linearform=f
xmass 1.0 1.0

#define preconditioner c -type=local -bilinearform=a -test #-block
define preconditioner c -type=direct -bilinearform=a -inverse=pardiso #-test
#define preconditioner c -type=bddc -bilinearform=a -test -block
#define preconditioner c -type=multigrid -bilinearform=a -test #-smoother=block

#numproc shapetester npst -gridfunction=u #-comp=0

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6

# for the evaluation of gradients
define bilinearform a_d -fespace=fescomp -nonassemble
laplace one -comp=1

numproc drawflux npdf -solution=u -bilinearform=a_d -label=grad

numproc visualization npviz -scalarfunction=u #-comp=0

# evaluate l2-error (difference between prescribed solution and numerical solution)
numproc xdifference npxd 
        -solution=u 
        -solution_n=solneg
        -solution_p=solpos
        -levelset=lset
        -intorder=3
